"""Metadata-only workflow for the compact O3 transition full-W sector leaf.

This external draft owns lifecycle metadata only.  It authenticates the saved
P2 parent and the static source-facts packet, freezes the input packet, and
defines the one future ``derive(source, diagnostic=None)`` entry point.  It
does not contain the W-sector mathematics, a provider, checker, reducer,
normalizer, or numerical engine.  Importing this file performs no filesystem,
Git, SymPy, kernel, or runner action.

The ``run`` action is deliberately written for a later controlled invocation;
this draft has not run it.  Its diagnostic callback writes each bare kernel
event immediately after the two manager prefix events, and its complete raw
result is persisted before validation.  Publication and resume operate on
those saved bytes and never rerun mathematics.
"""

from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import hashlib
import importlib.metadata as importlib_metadata
import importlib.util
import json
import os
from pathlib import Path
import subprocess
from typing import Any


REPO = Path(
    r"C:/Users/z5242/Documents/Codex/2026-09-13/"
    r"before-doing-any-project-work-verify/work/st073-transition-next"
)
WORK = Path(
    r"C:/Users/z5242/Documents/Codex/2026-10-16/"
    r"https-github-com-liang-techh-openai/work"
)

COMPONENT = "compact_O3_transition_W_sector_correlation"
INPUT_SCHEMA = "compact_O3_transition_W_sector_correlation_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_W_sector_correlation"
RESULT_SCHEMA_ALIASES = (
    RESULT_SCHEMA,
    "compact_O3_transition_W_sector_correlation_result_v1",
)
ENTERED_SCHEMA = "original_compact_O3_transition_W_sector_correlation_entered_v1"
DIAGNOSTIC_SCHEMA = "compact_O3_transition_W_sector_correlation_stage_diagnostic_v1"
PUBLICATION_SCHEMA = "original_compact_O3_transition_W_sector_correlation_publication_v1"
SEAL_SCHEMA = "compact_O3_transition_W_sector_correlation_final_seal_v1"

KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/"
    "original_compact_O3_transition_W_sector_correlation.py"
)
RUNNER_RELATIVE = Path(
    "scripts/run_original_compact_O3_transition_W_sector_correlation.py"
)
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_W_sector_correlation.json"
)
MILESTONE_RELATIVE = Path(
    "docs/CURRENT_ORIGINAL_COMPACT_O3_TRANSITION_W_SECTOR_CORRELATION_2026_10_27.md"
)
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
    MILESTONE_RELATIVE,
    *HANDOFF_RELATIVES,
)

CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = (
    SOURCE_CHECKPOINT,
    ENTERED_CHECKPOINT,
    RESULT_CHECKPOINT,
    DIAGNOSTIC_CHECKPOINT,
)

INPUT_SPEC = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"
RESUME_RECEIPT = WORK / f"{COMPONENT}_resume_receipt_2026_10_27.json"
SEAL_RECEIPT = WORK / f"{COMPONENT}_final_seal_receipt_2026_10_27.json"

PARENT_REPORT_DEFAULT = (
    REPO
    / "experiments/root_st073/"
    / "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_P2_pressure_gap.json"
)
PARENT_REPORT_SHA256 = (
    "2956baf8c82250ffd789714ceb1ee80170bd69e57f3c84d1d6b3b8889a1d1466"
)
PARENT_RESULT_DEFAULT = (
    REPO / "tmp/compact_O3_transition_P2_pressure_gap/result.json"
)
PARENT_RESULT_SHA256 = (
    "958ba5154ea7a532f33d2e48cafeb12209bf94e376b9a0cc5dafae83bfb676a7"
)
PARENT_COMMIT = "fec736e2a0ae6b55a59f7ab22f0e2037d7a5e4f3"
PARENT_PARENT_COMMIT = "b7687e8f187ab5e54552d9a3d28201d69a38476c"

SOURCE_FACTS_DEFAULT = WORK / "wsector_source_trace_2026_10_27_v2.json"
SOURCE_FACTS_SHA256 = (
    "eb0af8bce8ab9af4bb72f8bc8a3eb5e03f5a2a6db9688ab30f37104084003dcf"
)
SOURCE_FACTS_V1_SHA256 = (
    "7d11e08d8b52ef2eab5ff88b7b49e8087cd296b47c59a77d933d662b36475e56"
)
DESIGN_REVIEW_DEFAULT = WORK / "wsector_design_review_2026_10_27_v2.json"
DESIGN_REVIEW_SHA256 = (
    "a7797c06cebf04176d3c4cc5e714e87d9f0669c9f278fd3043cd951c24c90d4c"
)
PARAMETER_ALIAS_TRACE_DEFAULT = WORK / "wsector_parameter_alias_trace_2026_10_27.json"
PARAMETER_ALIAS_TRACE_SHA256 = (
    "b483e1bd0657fc2ad223f4153044875191b5beae7e2028e99f7a08b9d05ba112"
)
PARAMETER_ALIAS_TRACE_BYTES = 23979
PARENT_REVIEW_DEFAULT = (
    WORK / "review_compact_O3_transition_P2_pressure_gap_saved_2026_10_27_v2_receipt.json"
)
PARENT_REVIEW_SHA256 = (
    "31dfc0eb524bcc767a4d86766c719aeddc840497b1b9b32f50a64b2855ad876d"
)

PARENT_SOURCE_PIN_COUNT = 2234
PARENT_CHECKPOINT_COUNT = 252
CURRENT_CHECKPOINT_COUNT = 256
PARENT_FAILED_COUNT = 176
FAILED_COUNT_IF_NEW_FAILURE = PARENT_FAILED_COUNT + 7
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26

SYMPY_MODULE_DEFAULTS = (
    "sympy/__init__.py",
    "sympy/core/function.py",
    "sympy/core/expr.py",
)
SYMPY_EXPECTED_SHA256 = {
    "sympy/__init__.py": "4e9476348ba105feab28d82f5bcf6cdba2e3e84de6e059bbfe7a13728c0a4ab0",
    "sympy/core/function.py": "cb4bc47576940fcb613f75dab1956492882a5ce181d1cb362677c711ae033634",
    "sympy/core/expr.py": "1d1ff1ce1eb6ef714199d7cd6e8806ce45a8a9a23f93f09e20e9e99b2cf00ab4",
    "sympy-1.14.0.dist-info/RECORD": "36253d75a569bb0b282177139eea1e59fe6d6df0b5e364f570ca41d71473181b",
    "sympy-1.14.0.dist-info/METADATA": "b756c2fbfd5be05ac5bdb0ebca61f55618f30f633ed92d88a5687429313a7595",
}

EXPECTED_COUNTS = {
    "new_mathematical_proofs": 1,
    "new_algebra_family_calls": 1,
    "source_bound_sides": 2,
    "sector_families": 8,
    "numerator_jet_rows": 32,
    "cross_numerator_identity_orders": 4,
    "whole_W_gap_derivative_views": 4,
    "symbolic_diff_api_calls": 36,
    "symbolic_expand_api_calls": 4,
    "numeric_calls": 0,
    "source_provider_calls": 0,
    "old_math_calls": 0,
}
COUNT_KEYS = tuple(EXPECTED_COUNTS)
NEW_FLAGS = (
    "O3_transition_full_W_signed_sector_functions_defined",
    "O3_transition_full_W_cross_numerator_correlation_ordinary_Z3_proved",
)

# These are exact field names, not substring scans.  A descriptive Boolean
# containing words such as recurrence or provider is never treated as a call
# counter.
NONCLAIM_FLAGS = (
    "full_goal_complete",
    "physical_layer_completed",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "global_closure_proved",
    "time_order_recursion_claimed",
    "O3_transition_full_W_zero_claimed",
    "O3_transition_phase_equality_proved",
    "O3_transition_signed_gap_zero_claimed",
    "pointwise_parameter_identity_proved",
    "radius_identity_proved",
    "cap_error_reduction_claimed",
    "error_cap_reduced",
    "old_error_terms_removed",
    "old_errors_removed",
    "old_errors_removed_global",
    "global",
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
    "global_five_moment_closure",
    "numeric_field_values_present",
)

AUDIT_ZERO_KEYS = (
    "numeric_calls",
    "source_provider_calls",
    "old_math_calls",
    "provider_calls",
    "checker_calls",
    "normalizer_calls",
    "reducer_calls",
    "quadrature_calls",
    "recurrence_calls",
    "field_evaluations",
    "integral_evaluations",
)
ANCESTOR_ZERO_KEYS = (
    "ancestor_math_helper_calls",
    "ancestor_provider_calls",
    "ancestor_checker_calls",
    "ancestor_normalizer_calls",
    "ancestor_reducer_calls",
    "ancestor_quadrature_calls",
    "ancestor_recurrence_calls",
    "ancestor_numeric_calls",
)

RETAINED_MAP_NAMES = (
    "preserved_error_cap_metadata",
    "retained_construction_data",
    "retained_integral_term_inventory",
    "retained_cap_provenance",
)

KERNEL_CONTRACT = {
    "derive": "derive(source, diagnostic=None)",
    "result_schema": RESULT_SCHEMA,
    "kernel_result_metadata_only": False,
    "kernel_result_source_only": True,
    "generic_engine": "SymPy 1.14.0",
    "callback_prefix_events": 2,
    "callback_event_count": 41,
    "row_callback_event_count": 40,
    "algebra_entry_callback_count": 1,
    "sector_jet_callback_count": 32,
    "cross_identity_callback_count": 4,
    "whole_W_gap_views_are_unevaluated": True,
    "whole_W_gap_definition": "D_Z^j[(sum_of_8_cross_numerator_sectors)/(L_n*L_c)]",
    "held_fixed_chart": "R/S/delta held fixed in ordinary Z; independent across sides",
    "physical_common_point_identity": False,
    "numeric_or_physical_claims": False,
}

NORMALIZED_FORMULAS = {
    "A": "-V+(1-delta)*(m-Z*m_Z)/2",
    "B": "V*((1-delta)*Z*m+(1-Z**2)*m_Z)+2*delta*Z*e-(1-Z**2)*e_Z+2*(1+delta)*Z*(P0+p)-(1-Z**2)*(P0_Z+p_Z)",
    "L": "1-delta*Z**2",
    "W": "R*(A+S*B)/L",
}
KNOWN_AST_DIGESTS = (
    "e63ae3283447c56ae542f03f2ab1480bb7e96d1624cc1b83a3cb53b88354f11c",
    "fbfe367ab432bd721d16129adf80595e0f512fe7b8a8bff62cb27378a29bb005",
    "ab33333c0bff2bd25ce8b8a0aa0472f656d43d1628b350103c620311d8205aaa",
    "3ebe8c6fa524fb0c19e0982313fa02c1957b6f919cd9b4ea793988c51560d1fe",
    "ae7ea5b8440f9f654dc556acfa864a7df15cacb594642470eae6efb4b7e83da3",
    "58101254f7a8609a45ed3ee832c47e7bd185338c326795f334e49b3635422648",
    "984da1fa05f96f9d1561df51e67b26ef77f0cf223806d96c8ededadb9d08b87b",
    "f68ccb8f07e6244113f92404aa30cf904176df03bda8f6508600182c4efe0dfc",
)

REQUIRED_SOURCE_FUNCTION_KEYS = (
    "native_generic_functions",
    "native_generic_roots_handoff",
    "native_outer_call_route",
    "native_outer_transport_geometry",
    "native_quotient_cell_port",
    "rebuilt_C3_recover_inputs",
    "rebuilt_C3_source_quotients_port",
    "rebuilt_compile_recovery_ast_anchor",
)


def require(condition: Any, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_path(path: Path) -> str:
    return sha_bytes(Path(path).read_bytes())


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


def clone(value: Any) -> Any:
    return json.loads(encode(value).decode("utf-8"))


def brief_safe(value: Any) -> Any:
    """Compact BRIEF aliases without changing publication/source artifacts."""
    if isinstance(value, dict):
        return {
            key: brief_safe(child)
            for key, child in value.items()
            if key != "raw_bytes_base64"
        }
    if isinstance(value, list):
        return [brief_safe(child) for child in value]
    return value


def wire_path(path: Path | str) -> str:
    return Path(path).as_posix()


def portable(value: Any) -> Any:
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            portable_key = key.replace("\\", "/") if isinstance(key, str) else key
            require(portable_key not in result, "Portable metadata key collision")
            result[portable_key] = portable(child)
        return result
    if isinstance(value, list):
        return [portable(child) for child in value]
    if isinstance(value, str) and ("\\" in value or (len(value) > 2 and value[1] == ":")):
        return value.replace("\\", "/")
    return value


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


def valid_git_revision(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 40:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


def resolve_path(value: str | Path, base: Path = WORK) -> Path:
    """Resolve absolute and relative pins, including root_st073 source pins."""
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    candidates = (
        base / candidate,
        REPO / candidate,
        REPO / "experiments/root_st073" / candidate,
        WORK / candidate,
    )
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
            old = result.get(path)
            require(old is None or old == digest, f"Hash collision changed: {path}")
            result[path] = digest
    return result


def artifact_entry(
    path: Path,
    digest: str | None = None,
    embed: bool = False,
    outside_science_pins: bool = False,
) -> dict:
    path = Path(path)
    raw = read_bytes(path, digest)
    entry = {
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


def bounded_parent_view(report: dict) -> dict:
    """Keep headers/counts/flags only; never copy nested full reports or base64."""
    keys = (
        "schema",
        "component",
        "attempt_status",
        "passed",
        "metadata_only",
        "source_only",
        "parent_commit",
        "parent_head",
        "source_pin_count",
        "checkpoint_count",
        "preserved_failed_artifact_count",
        "parent_protected_file_count",
        "parent_dirty_path_count",
        "expected_parent_false_flag_count",
        "actual_counts",
        "new_flags",
        "full_goal_complete",
        "parent_false_flags",
        "parent_false_flag_values",
    )
    view = portable({key: report[key] for key in keys if key in report})
    view["parent_report_is_authenticated_view"] = True
    view["complete_parent_report_raw_sha256"] = PARENT_REPORT_SHA256
    view["omitted_fields_saved_in_complete_parent_report"] = sorted(
        set(report).difference(keys)
    )
    require(len(encode(view)) < 100_000, "Bounded parent view exceeds 100 KB")
    return view


def load_parent(args: argparse.Namespace) -> dict:
    report_path = resolve_path(args.parent_report or PARENT_REPORT_DEFAULT, WORK)
    report_sha = (args.parent_report_sha256 or PARENT_REPORT_SHA256).lower()
    require(report_sha == PARENT_REPORT_SHA256, "Parent report SHA is not authenticated")
    report_raw, report = read_json(report_path, report_sha)
    require(
        report.get("schema") == "original_compact_O3_transition_P2_pressure_gap_publication_v1",
        "Wrong parent report schema",
    )
    require(
        report.get("component") == "compact_O3_transition_P2_pressure_gap",
        "Wrong parent report component",
    )
    require(report.get("passed") is True, "Parent publication did not pass")
    require(report.get("metadata_only") is True, "Parent publication is not metadata-only")
    require(report.get("source_only") is True, "Parent publication is not source-only")
    require(report.get("parent_commit") == PARENT_PARENT_COMMIT, "Parent report anchor changed")
    require(report.get("parent_head") == PARENT_PARENT_COMMIT, "Parent report head changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Parent source-pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Parent checkpoint count changed")
    require(
        report.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT,
        "Parent failed-artifact count changed",
    )
    source_pins = auth_map(report.get("source_pins"), "Parent source pins", PARENT_SOURCE_PIN_COUNT)
    parent_checkpoints = auth_map(
        report.get("parent_checkpoint_hashes"),
        "Parent checkpoint hashes",
        PARENT_CHECKPOINT_COUNT - 4,
    )
    current_checkpoints = auth_map(report.get("checkpoint_hashes"), "Parent current checkpoint hashes", 4)
    for name, digest in merge_hashes(parent_checkpoints, current_checkpoints).items():
        read_bytes(resolve_path(name, REPO), digest)
    protected = auth_map(
        report.get("protected_file_hashes"), "Parent protected files", PARENT_PROTECTED_COUNT
    )
    dirty = auth_map(
        report.get("legacy_dirty_path_hashes"), "Parent dirty paths", PARENT_DIRTY_COUNT
    )
    failed = auth_map(
        report.get("preserved_failed_artifact_hashes_union"),
        "Parent failed artifacts",
        PARENT_FAILED_COUNT,
    )
    parent_flags = report.get("parent_false_flags")
    parent_values = report.get("parent_false_flag_values")
    require(isinstance(parent_flags, list), "Parent false-flag ledger is not a list")
    require(len(parent_flags) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag ledger count changed")
    require(len(set(parent_flags)) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag ledger duplicated")
    require(isinstance(parent_values, dict), "Parent false-flag values are missing")
    require(parent_values == {key: False for key in parent_flags}, "Parent false values changed")
    for name in RETAINED_MAP_NAMES:
        require(isinstance(report.get(name), dict), f"Full parent retained map missing: {name}")

    result_path = resolve_path(PARENT_RESULT_DEFAULT, REPO)
    result_raw, result = read_json(result_path, PARENT_RESULT_SHA256)
    require(report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256, "Parent result SHA changed")
    require(report.get("full_result") == result, "Parent full result is not exact")
    require(result.get("schema") == "original_compact_O3_transition_P2_pressure_gap", "Parent result schema changed")
    require(result.get("component") == "compact_O3_transition_P2_pressure_gap", "Parent result component changed")
    require(result.get("passed") is True, "Saved parent result did not pass")
    return {
        "report": report,
        "report_raw": report_raw,
        "report_meta": {
            "path": wire_path(report_path),
            "sha256": PARENT_REPORT_SHA256,
            "bytes": len(report_raw),
            "schema": report.get("schema"),
            "component": report.get("component"),
        },
        "report_view": bounded_parent_view(report),
        "result": result,
        "result_raw": result_raw,
        "result_meta": {
            "path": wire_path(result_path),
            "sha256": PARENT_RESULT_SHA256,
            "bytes": len(result_raw),
            "schema": result.get("schema"),
            "component": result.get("component"),
        },
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


def load_source_facts(path: Path, expected_sha: str) -> tuple[bytes, dict]:
    raw, facts = read_json(path, expected_sha)
    require(facts.get("packet_type") == "wsector_source_trace_v2", "Wrong W-sector source packet")
    require(facts.get("status") == "static_source_ast_trace_only", "Source packet is not static-only")
    require(facts.get("v1_artifact", {}).get("sha256") == SOURCE_FACTS_V1_SHA256, "V1 source trace SHA changed")
    raw_inputs = facts.get("raw_inputs")
    require(isinstance(raw_inputs, dict), "Source packet raw-input metadata is missing")
    require(raw_inputs.get("parent_report", {}).get("sha256") == PARENT_REPORT_SHA256, "Source parent report pin changed")
    require(raw_inputs.get("saved_result", {}).get("sha256") == PARENT_RESULT_SHA256, "Source saved-result pin changed")
    scope = facts.get("scope", {})
    v2_scope = facts.get("v2_scope", {})
    require(scope.get("source_bytes_and_AST_only") is True, "Source packet is not bytes/AST-only")
    require(scope.get("project_modules_executed") is False, "Source packet executed project modules")
    require(scope.get("project_modules_imported") is False, "Source packet imported project modules")
    require(v2_scope.get("new_source_consumer_is_metadata_only") is True, "Source consumer is not metadata-only")
    require(v2_scope.get("project_modules_executed") is False, "V2 source packet executed project modules")
    require(v2_scope.get("project_modules_imported") is False, "V2 source packet imported project modules")
    require(isinstance(facts.get("v2_normalized_formula_facts"), dict), "Normalized source facts are missing")
    anchors = facts.get("v2_source_anchors")
    require(isinstance(anchors, dict), "Source AST anchors are missing")
    coordinate_binding = anchors.get("native_added_transport_and_coordinate_bindings", {}).get("coordinate_Z_binding", {})
    require(coordinate_binding.get("targets") == ["z"], "Source ordinary-Z binding target changed")
    source_functions = facts.get("source_functions")
    require(isinstance(source_functions, dict), "Source function IR is missing")
    require(
        set(REQUIRED_SOURCE_FUNCTION_KEYS).issubset(source_functions),
        "Source function IR lacks an actual-port or AST anchor",
    )
    signed_sectors = facts.get("full_signed_axial_sectors")
    require(isinstance(signed_sectors, dict), "Full signed sector facts are missing")
    require(all(isinstance(signed_sectors.get(side), dict) for side in ("native", "rebuilt_C3")), "Both source sides are required")
    values = []

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for key, child in node.items():
                values.append((str(key), child))
                walk(child)
        elif isinstance(node, list):
            for child in node:
                walk(child)

    walk(facts)
    text_values = {child for _key, child in values if isinstance(child, str)}
    for name, formula in NORMALIZED_FORMULAS.items():
        require(formula in text_values or f"{name}={formula}" in text_values, f"Normalized formula pin missing: {name}")
    digest_values = {child for _key, child in values if isinstance(child, str) and valid_sha(child)}
    require(set(KNOWN_AST_DIGESTS).issubset(digest_values), "Source AST digest set is incomplete")
    function_values = {child for key, child in values if key == "function" and isinstance(child, str)}
    require({"generic_functions", "recover_inputs"}.issubset(function_values), "Typed producer function anchors are incomplete")
    return raw, facts


def load_design_review() -> tuple[bytes, dict]:
    raw, review = read_json(DESIGN_REVIEW_DEFAULT, DESIGN_REVIEW_SHA256)
    require(review.get("schema") == "wsector_design_review.v2_addendum", "Wrong W-sector design review schema")
    require(review.get("status") == "static_advisory_only_no_execution", "Design review is not static-only")
    require(review.get("receipt_state", {}).get("kernel_executed") is False, "Design review records kernel execution")
    require(review.get("receipt_state", {}).get("proof_executed") is False, "Design review records proof execution")
    require(review.get("chart_vs_physical_scope", {}).get("non_identity_statement"), "Design review physical nonidentity guard missing")
    chart = review.get("held_fixed_chart_contract", {})
    require(
        chart.get("held_fixed_scalars") == ["Rn", "Rc", "Sn", "Sc", "dn", "dc"],
        "Held-fixed scalar ports changed",
    )
    require(review.get("jet_exports_and_planned_counts", {}).get("planned_derivative_calls") == 36, "Design review diff count changed")
    require(review.get("jet_exports_and_planned_counts", {}).get("planned_expand_calls") == 4, "Design review expand count changed")
    return raw, review


def parameter_alias_trace_artifact() -> dict:
    """Pin static alias evidence outside the mathematical input packet."""
    raw, trace = read_json(PARAMETER_ALIAS_TRACE_DEFAULT, PARAMETER_ALIAS_TRACE_SHA256)
    require(len(raw) == PARAMETER_ALIAS_TRACE_BYTES, "Parameter alias trace byte count changed")
    require(trace.get("schema") == "wsector_parameter_alias_trace_static_v1", "Wrong parameter alias trace schema")
    require(trace.get("runtime_model_metadata") is None, "Parameter alias trace has runtime metadata")
    audit = trace.get("static_audit", {})
    require(isinstance(audit, dict), "Parameter alias static audit is missing")
    require(audit.get("project_modules_imported") is False, "Parameter alias trace imported project modules")
    require(audit.get("project_modules_called") is False, "Parameter alias trace called project modules")
    require(audit.get("math_or_helper_consumer_calls") == 0, "Parameter alias trace contains math/helper calls")
    status = trace.get("alias_and_theorem_status", {})
    require(isinstance(status, dict), "Parameter alias theorem status is missing")
    for field in (
        "native_to_rebuilt_same_source_object_Pstar",
        "native_to_rebuilt_same_source_object_delta",
        "native_to_rebuilt_same_source_object_radius",
        "explicit_saved_alias_theorem_for_delta_S",
        "explicit_saved_alias_theorem_for_delta_delta",
    ):
        require(status.get(field) is False, f"Parameter alias status changed: {field}")
    require(trace.get("next_action_readiness", {}).get("no_math_identity_claim") is True, "Parameter alias trace makes a math identity claim")
    return artifact_entry(
        PARAMETER_ALIAS_TRACE_DEFAULT,
        PARAMETER_ALIAS_TRACE_SHA256,
        embed=True,
        outside_science_pins=True,
    )


def sympy_raw_pins(module_files: list[str]) -> dict:
    """Pin raw distribution metadata and exactly three files without importing SymPy."""
    require(len(module_files) == 3, "Exactly three SymPy module files are required")
    distribution = importlib_metadata.distribution("sympy")
    require(distribution.version == "1.14.0", "SymPy version is not the frozen 1.14.0")
    dist_files = list(distribution.files or ())
    metadata_candidates = [
        distribution.locate_file(item)
        for item in dist_files
        if str(item).replace("\\", "/").endswith(".dist-info/METADATA")
    ]
    record_candidates = [
        distribution.locate_file(item)
        for item in dist_files
        if str(item).replace("\\", "/").endswith(".dist-info/RECORD")
    ]
    require(len(metadata_candidates) == 1 and len(record_candidates) == 1, "SymPy METADATA/RECORD pin is ambiguous")
    metadata_path = metadata_candidates[0].resolve()
    record_path = record_candidates[0].resolve()
    dist_root = Path(distribution.locate_file(".")).resolve()
    metadata_key = metadata_path.relative_to(dist_root).as_posix()
    record_key = record_path.relative_to(dist_root).as_posix()
    require(metadata_key in SYMPY_EXPECTED_SHA256, "Unexpected SymPy METADATA path")
    require(record_key in SYMPY_EXPECTED_SHA256, "Unexpected SymPy RECORD path")
    entries = [
        artifact_entry(
            metadata_path,
            SYMPY_EXPECTED_SHA256[metadata_key],
            embed=True,
            outside_science_pins=True,
        ),
        artifact_entry(
            record_path,
            SYMPY_EXPECTED_SHA256[record_key],
            embed=True,
            outside_science_pins=True,
        ),
    ]
    seen = set()
    for module_file in module_files:
        path = Path(module_file)
        path = path if path.is_absolute() else distribution.locate_file(path)
        path = path.resolve()
        require(path.is_file(), f"SymPy module file is missing: {module_file}")
        require(path not in seen, "SymPy module file repeated")
        seen.add(path)
        module_key = path.relative_to(dist_root).as_posix()
        require(module_key in SYMPY_EXPECTED_SHA256, f"Unexpected SymPy module file: {module_key}")
        entries.append(
            artifact_entry(
                path,
                SYMPY_EXPECTED_SHA256[module_key],
                embed=True,
                outside_science_pins=True,
            )
        )
    return {"version": distribution.version, "artifacts": entries}


def build_input_spec(
    parent: dict,
    facts_path: Path,
    facts: dict,
    facts_raw: bytes,
    design_review: dict,
    design_review_raw: bytes,
    review: dict,
    sympy: dict,
) -> tuple[bytes, dict]:
    value = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_report": parent["report_meta"],
        "parent_result_metadata": parent["result_meta"],
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_result": clone(parent["result"]),
        "sector_source_facts": clone(facts),
        "sector_source_facts_path": wire_path(facts_path),
        "sector_source_facts_raw_sha256": sha_bytes(facts_raw),
        "sector_source_facts_artifact": artifact_entry(
            facts_path, sha_bytes(facts_raw), embed=True, outside_science_pins=True
        ),
        "design_review_artifact": artifact_entry(
            DESIGN_REVIEW_DEFAULT,
            sha_bytes(design_review_raw),
            embed=True,
            outside_science_pins=True,
        ),
        "design_review": clone(design_review),
        "review": review,
        "sympy_raw_pins": clone(sympy),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": list(NEW_FLAGS),
        "kernel_contract": clone(KERNEL_CONTRACT),
        "derive_contract": KERNEL_CONTRACT["derive"],
        "expected_kernel_result_metadata_only": KERNEL_CONTRACT["kernel_result_metadata_only"],
        "whole_W_gap_definition": KERNEL_CONTRACT["whole_W_gap_definition"],
        "whole_W_gap_views_unevaluated": KERNEL_CONTRACT["whole_W_gap_views_are_unevaluated"],
        "conditional_source_functions": True,
        "physical_layer_completed": False,
        "science_reexecuted": False,
        "metadata_only_assembly": True,
    }
    raw = write_new(INPUT_SPEC, value)
    value["raw_sha256"] = sha_bytes(raw)
    return raw, value


def build_source(
    parent: dict,
    facts_path: Path,
    input_raw: bytes,
    facts: dict,
    facts_raw: bytes,
    design_review: dict,
    design_review_raw: bytes,
    review: dict,
    sympy: dict,
    kernel_path: Path,
    runner_path: Path,
) -> dict:
    kernel_raw = read_bytes(kernel_path)
    runner_raw = read_bytes(runner_path)
    v1_path = resolve_path(facts.get("v1_artifact", {}).get("path", ""), WORK)
    v1_raw = read_bytes(v1_path, SOURCE_FACTS_V1_SHA256)
    additions = {
        parent["report_meta"]["path"]: PARENT_REPORT_SHA256,
        parent["result_meta"]["path"]: PARENT_RESULT_SHA256,
        wire_path(INPUT_SPEC): sha_bytes(input_raw),
        wire_path(facts_path): sha_bytes(facts_raw),
        wire_path(v1_path): SOURCE_FACTS_V1_SHA256,
        wire_path(DESIGN_REVIEW_DEFAULT): sha_bytes(design_review_raw),
        wire_path(kernel_path): sha_bytes(kernel_raw),
        wire_path(runner_path): sha_bytes(runner_raw),
    }
    for entry in sympy["artifacts"]:
        additions[entry["path"]] = entry["sha256"]
    source_union = merge_hashes(parent["source_pins"], additions)
    source = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_parent_commit": PARENT_PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "parent_report_view": parent["report_view"],
        "parent_result": clone(parent["result"]),
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_result_data_only": True,
        "parent_passed": True,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "parent_checkpoint_hashes": merge_hashes(parent["parent_checkpoints"], parent["current_checkpoints"]),
        "parent_source_pins": parent["source_pins"],
        "parent_protected_hashes": parent["protected"],
        "parent_dirty_hashes": parent["dirty"],
        "parent_failed_hashes": parent["failed"],
        "parent_false_flags": parent["false_flags"],
        "parent_false_flag_values": parent["false_flag_values"],
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "preserved_error_cap_metadata": parent["retained"]["preserved_error_cap_metadata"],
        "retained_construction_data": parent["retained"]["retained_construction_data"],
        "retained_integral_term_inventory": parent["retained"]["retained_integral_term_inventory"],
        "retained_cap_provenance": parent["retained"]["retained_cap_provenance"],
        "sector_source_facts": clone(facts),
        "sector_source_facts_raw_sha256": sha_bytes(facts_raw),
        "sector_source_facts_artifact": artifact_entry(
            facts_path, sha_bytes(facts_raw), outside_science_pins=True
        ),
        "sector_source_facts_v1_artifact": artifact_entry(
            v1_path, SOURCE_FACTS_V1_SHA256, outside_science_pins=True
        ),
        "design_review": clone(design_review),
        "design_review_artifact": artifact_entry(
            DESIGN_REVIEW_DEFAULT,
            sha_bytes(design_review_raw),
            outside_science_pins=True,
        ),
        "review": review,
        "sympy_raw_pins": clone(sympy),
        "input_spec_metadata": {
            "path": wire_path(INPUT_SPEC),
            "sha256": sha_bytes(input_raw),
            "bytes": len(input_raw),
            "raw_readonly": True,
        },
        "input_raw_sha256": sha_bytes(input_raw),
        "kernel_path": wire_path(kernel_path),
        "kernel_sha256": sha_bytes(kernel_raw),
        "runner_path": wire_path(runner_path),
        "runner_sha256": sha_bytes(runner_raw),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": list(NEW_FLAGS),
        "kernel_contract": clone(KERNEL_CONTRACT),
        "expected_kernel_result_metadata_only": KERNEL_CONTRACT["kernel_result_metadata_only"],
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "raw_sha_pins": source_union,
        "source_packet_defined": True,
        "source_packet_definition_accepted": True,
        "source_only": True,
        "metadata_only": True,
        "conditional_source_functions": True,
        "physical_layer_completed": False,
        "pointwise_identity_proved": False,
        "cap_error_reduction_claimed": False,
        "accepted_progress_physical_baseline_preserved": True,
        "runtime_model_metadata": None,
        "math_reexecuted": False,
        "old_math_reexecuted": False,
        "no_runner_execution": True,
    }
    return source


def load_source() -> tuple[bytes, dict]:
    raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == INPUT_SCHEMA, "Source packet schema changed")
    require(source.get("component") == COMPONENT, "Source packet component changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source parent commit changed")
    require(source.get("parent_report_raw_sha256") == PARENT_REPORT_SHA256, "Source parent report SHA changed")
    require(source.get("parent_result_raw_sha256") == PARENT_RESULT_SHA256, "Source parent result SHA changed")
    require(source.get("sector_source_facts_raw_sha256") == SOURCE_FACTS_SHA256, "Sector source-facts SHA changed")
    require(source.get("sector_source_facts_path"), "Sector source-facts path is missing")
    require(source.get("design_review_artifact", {}).get("sha256") == DESIGN_REVIEW_SHA256, "Design review SHA changed")
    require(source.get("kernel_contract") == KERNEL_CONTRACT, "Kernel contract changed")
    require(source.get("parent_false_flags") and len(source["parent_false_flags"]) == PARENT_FALSE_FLAG_COUNT, "Parent ledger changed")
    require(source.get("parent_false_flag_values") == {key: False for key in source["parent_false_flags"]}, "Parent ledger values changed")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    for name, digest in source.get("source_union_hashes", {}).items():
        read_bytes(resolve_path(name, WORK), digest)
    for name in RETAINED_MAP_NAMES:
        require(isinstance(source.get(name), dict), f"Retained parent map missing: {name}")
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


def diagnostic_summary() -> dict:
    lines = DIAGNOSTIC_CHECKPOINT.read_bytes().splitlines()
    events = []
    for line in lines:
        require(line, "Empty diagnostic line")
        event = json.loads(line.decode("utf-8"))
        require(isinstance(event, dict), "Diagnostic line is not an object")
        events.append(event)
    if events:
        require(events[0].get("event") == "manager_entry_started", "Diagnostic prefix one changed")
    if len(events) > 1:
        require(
            events[1].get("event") in ("manager_kernel_imported", "manager_kernel_or_derive_failed"),
            "Diagnostic prefix two changed",
        )
    callback_events = events[2:]
    event_counts = {}
    for event in callback_events:
        name = event.get("event")
        event_counts[name] = event_counts.get(name, 0) + 1
    return {
        "event_count": len(events),
        "manager_prefix_event_count": min(len(events), 2),
        "manager_kernel_imported": len(events) > 1 and events[1].get("event") == "manager_kernel_imported",
        "callback_event_count": len(callback_events),
        "row_callback_event_count": sum(
            event_counts.get(name, 0)
            for name in ("numerator_jet_row", "direct_source_jet_row", "cross_numerator_identity_row")
        ),
        "callback_event_counts": event_counts,
        "callback_events_are_bare": all("manager_event" not in event for event in callback_events),
        "callback_event_raw_rows": [clone(event) for event in callback_events],
    }


def import_kernel(path: Path):
    spec = importlib.util.spec_from_file_location("_wsector_kernel", path)
    require(spec is not None and spec.loader is not None, "Kernel module spec unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "Kernel derive(source, diagnostic=None) is missing")
    return module


def failure_result(
    error: BaseException,
    derive_entry_calls: int,
    derive_calls: int,
    kernel_import_calls: int,
) -> dict:
    partial_counts = {name: 0 for name in COUNT_KEYS}
    partial_flags = {name: False for name in NEW_FLAGS}
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "failed",
        "passed": False,
        "metadata_only": False,
        "source_only": True,
        "conditional_on_source_port_semantics": True,
        "physical_layer_completed": False,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_import_calls == 1,
        "error_type": type(error).__name__,
        "error_message": str(error),
        "actual_counts": partial_counts,
        "actual_api_counts": dict(partial_counts),
        "new_flags": partial_flags,
        "flags": dict(partial_flags),
        "runtime_model_metadata": None,
    }


def validate_result(
    result: dict,
    source: dict,
    derive_entry_calls: int | None = None,
    derive_calls: int | None = None,
    diagnostic_callback_count: int | None = None,
    diagnostic_manager_kernel_imported: bool | None = None,
    diagnostic_event_counts: dict | None = None,
    diagnostic_row_callback_count: int | None = None,
) -> dict:
    violations = []
    if not isinstance(result, dict):
        violations.append("result is not an object")
        return {"passed": False, "protocol_violations": violations, "actual_counts": None, "new_flags": {}}
    result_passed = result.get("passed", result.get("status") == "ok")
    if not isinstance(result_passed, bool):
        violations.append("passed is not an explicit Boolean")
    if result.get("schema") not in RESULT_SCHEMA_ALIASES:
        violations.append("result schema changed")
    if result.get("component") != COMPONENT:
        violations.append("result component changed")
    if result.get("metadata_only") is not False:
        violations.append("kernel metadata_only marker is not false")
    if result.get("source_only") is not True:
        violations.append("source_only is not true")
    if result.get("physical_layer_completed") is True:
        violations.append("physical layer was claimed")
    component_result = result.get("component_result")
    if component_result is not None:
        if not isinstance(component_result, dict):
            violations.append("component_result is not an object")
        else:
            if result.get("component_result_exact") is not True:
                violations.append("component_result exact-preservation marker is missing")
            if component_result.get("schema") not in RESULT_SCHEMA_ALIASES:
                violations.append("component_result schema changed")
            if component_result.get("component") != COMPONENT:
                violations.append("component_result component changed")
            if component_result.get("metadata_only") is not False:
                violations.append("component_result metadata_only marker is not false")
            if component_result.get("source_only") is not True:
                violations.append("component_result source_only marker is not true")
            if result_passed and component_result.get("passed") is not True:
                violations.append("component_result passed marker is not true")
    conditional = result.get("conditional_source_functions", result.get("conditional_on_source_port_semantics"))
    if conditional is not True:
        violations.append("conditional source-function marker missing")
    if result_passed:
        formal_ports = result.get("formal_ports")
        if not isinstance(formal_ports, dict):
            violations.append("successful result lacks formal_ports")
        else:
            common_ports = formal_ports.get("common")
            if not isinstance(common_ports, dict) or common_ports.get("Z") != "Z":
                violations.append("successful result ordinary-Z machine coordinate changed")
            if formal_ports.get("parameters_held_fixed") != ["Rn", "Rc", "Sn", "Sc", "dn", "dc"]:
                violations.append("successful result held-fixed scalar ports changed")
    counts = result.get("actual_counts", result.get("actual_api_counts"))
    if result_passed:
        if counts != EXPECTED_COUNTS:
            violations.append("actual_counts do not exactly match the twelve-key contract")
    elif not isinstance(counts, dict):
        violations.append("failed result lacks partial actual_counts")
    else:
        unknown_counts = set(counts).difference(COUNT_KEYS)
        if unknown_counts:
            violations.append("failed result contains unknown count keys")
        for name, value in counts.items():
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                violations.append(f"failed result count is not a nonnegative integer: {name}")
    flags = result.get("new_flags", result.get("flags"))
    if result_passed:
        if not isinstance(flags, dict) or any(flags.get(name) is not True for name in NEW_FLAGS):
            violations.append("new W-sector flags are not both true")
    elif not isinstance(flags, dict) or any(not isinstance(value, bool) for value in flags.values()):
        violations.append("failed result flags are not Boolean metadata")
    for name in ("parent_false_flags", "inherited_false_flags"):
        if name in result and result[name] != source["parent_false_flags"]:
            violations.append(f"inherited false-flag list changed: {name}")
    for name in ("parent_false_flag_values", "inherited_false_flag_values"):
        if name in result and result[name] != source["parent_false_flag_values"]:
            violations.append(f"inherited false-flag values changed: {name}")
    for name in RETAINED_MAP_NAMES:
        if name in result and result[name] != source[name]:
            violations.append(f"retained parent map changed: {name}")
    for name in NONCLAIM_FLAGS:
        if name in result and result[name] is not False:
            violations.append(f"forbidden claim is true: {name}")
    for name in COUNT_KEYS:
        if isinstance(result.get(name), (int, float)):
            expected = EXPECTED_COUNTS[name] if result_passed else counts.get(name)
            if expected is not None and result[name] != expected:
                violations.append(f"top-level count disagrees: {name}")
    audit = result.get("audit")
    if isinstance(audit, dict):
        for name in AUDIT_ZERO_KEYS:
            if name in audit and audit[name] != 0:
                violations.append(f"audit count is nonzero: {name}")
    for name in ANCESTOR_ZERO_KEYS:
        if name in result and result[name] != 0:
            violations.append(f"ancestor count is nonzero: {name}")
    if "whole_W_gap_derivative_views" in result:
        views = result["whole_W_gap_derivative_views"]
        if not isinstance(views, list):
            violations.append("whole W gap views are not a list")
        elif result_passed and len(views) != 4:
            violations.append("whole W gap view count is not four")
        elif result_passed:
            for view in views:
                if not isinstance(view, dict):
                    violations.append("whole W gap view is not an object")
                    continue
                if view.get("source_conditional") is not True:
                    violations.append("whole W gap view lost source-conditional marker")
                if view.get("actual_symbolic_diff_calls") not in (None, 0):
                    violations.append("whole W gap view performed a symbolic diff")
                argument = str(view.get("argument", ""))
                if "Ln*Lc" not in argument:
                    violations.append("whole W gap view is not the whole quotient view")
    for name in ("whole_W_gap_views_unevaluated", "whole_W_gap_derivative_views_unevaluated"):
        if name in result and result[name] is not True:
            violations.append(f"whole W gap views were evaluated: {name}")
    if "whole_W_gap_definition" in result and result["whole_W_gap_definition"] != "D_Z^j[(sum_of_8_cross_numerator_sectors)/(L_n*L_c)]":
        violations.append("whole W gap definition changed")
    coordinate = result.get("kernel_workflow_coordinate")
    if coordinate is not None:
        if not isinstance(coordinate, dict):
            violations.append("kernel workflow coordinate is not an object")
        else:
            if coordinate.get("ordinary_Z_variable") != "Z":
                violations.append("ordinary-Z coordinate changed")
            if coordinate.get("held_fixed_scalars") != ["Rn", "Rc", "Sn", "Sc", "dn", "dc"]:
                violations.append("held-fixed scalar ports changed")
    claims = result.get("claims")
    if isinstance(claims, dict):
        for name in (
            "actual_numeric_jets",
            "source_W_gap_zero",
            "P2_gap_zero",
            "phase_equality",
            "radius_equality",
            "Pstar_equality",
            "delta_equality",
            "cap_reduction",
            "error_reduction",
            "physical_closure",
            "global_closure",
            "temporal_n_recursion",
        ):
            if name in claims and claims[name] is not False:
                violations.append(f"forbidden claim is true: claims.{name}")
    preserved = result.get("preserved_metadata")
    if isinstance(preserved, dict):
        for name in ("old_errors_removed", "old_error_terms_removed", "old_errors_removed_global", "global"):
            if name in preserved and preserved[name] is not False:
                violations.append(f"retained-error metadata changed: preserved_metadata.{name}")
        if "quiet_power4_retained" in preserved and preserved["quiet_power4_retained"] is not True:
            violations.append("quiet_power4 retained marker changed")
    if derive_entry_calls is not None:
        if result_passed and derive_entry_calls != 1:
            violations.append("successful derive_entry_calls is not one")
        if not result_passed and derive_entry_calls not in (0, 1):
            violations.append("failed derive_entry_calls is not a measured zero/one")
    if derive_calls is not None:
        if result_passed and derive_calls != 1:
            violations.append("successful derive_calls is not one")
        if not result_passed and derive_calls not in (0, 1):
            violations.append("failed derive_calls is not a measured zero/one")
    if "kernel_import_calls" in result:
        import_calls = result.get("kernel_import_calls")
        if isinstance(import_calls, bool) or not isinstance(import_calls, int) or import_calls not in (0, 1):
            violations.append("kernel_import_calls is not a measured zero/one")
        elif result_passed and import_calls != 1:
            violations.append("successful kernel_import_calls is not one")
    if "kernel_imported" in result:
        imported = result.get("kernel_imported")
        if not isinstance(imported, bool):
            violations.append("kernel_imported is not a Boolean")
        elif isinstance(result.get("kernel_import_calls"), int) and not isinstance(result.get("kernel_import_calls"), bool) and result.get("kernel_import_calls") in (0, 1) and imported is not (result.get("kernel_import_calls") == 1):
            violations.append("kernel_imported disagrees with kernel_import_calls")
    if diagnostic_callback_count is not None and result_passed is True and diagnostic_callback_count != 41:
        violations.append("diagnostic callback event count is not 41")
    if diagnostic_row_callback_count is not None and result_passed is True and diagnostic_row_callback_count != 40:
        violations.append("diagnostic row callback count is not 40")
    if diagnostic_manager_kernel_imported is not None and result_passed is True and not diagnostic_manager_kernel_imported:
        violations.append("successful result lacks kernel-import prefix")
    if diagnostic_event_counts is not None and result_passed is True:
        for event_name, expected in {
            "algebra_entry": 1,
            "numerator_jet_row": 32,
            "direct_source_jet_row": 4,
            "cross_numerator_identity_row": 4,
        }.items():
            if diagnostic_event_counts.get(event_name) != expected:
                violations.append(f"diagnostic event count changed: {event_name}")
    return {
        "passed": not violations and result_passed is True,
        "protocol_violations": violations,
        "actual_counts": clone(counts) if isinstance(counts, dict) else None,
        "new_flags": clone(flags) if isinstance(flags, dict) else {},
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "diagnostic_callback_count": diagnostic_callback_count,
        "diagnostic_row_callback_count": diagnostic_row_callback_count,
        "diagnostic_manager_kernel_imported": diagnostic_manager_kernel_imported,
        "diagnostic_event_counts": clone(diagnostic_event_counts) if isinstance(diagnostic_event_counts, dict) else None,
        "metadata_only": True,
        "source_only": True,
        "conditional_source_functions": conditional is True,
        "physical_layer_completed": False,
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": list(NEW_FLAGS),
    }


def prepare(args: argparse.Namespace) -> dict:
    require(not SOURCE_CHECKPOINT.exists(), "Source packet already exists; do not reenter")
    require(not INPUT_SPEC.exists(), "Input specification already exists; do not rewrite")
    require(all(not path.exists() for path in CHECKPOINTS), "A W-sector checkpoint already exists; do not reenter")
    parent = load_parent(args)
    facts_path = resolve_path(args.source_facts or SOURCE_FACTS_DEFAULT, WORK)
    facts_sha = (args.source_facts_sha256 or SOURCE_FACTS_SHA256).lower()
    require(facts_sha == SOURCE_FACTS_SHA256, "Source-facts SHA is not authenticated")
    facts_raw, facts = load_source_facts(facts_path, facts_sha)
    design_review_raw, design_review = load_design_review()
    review_path = resolve_path(args.review or PARENT_REVIEW_DEFAULT, WORK)
    review_sha = (args.review_sha256 or PARENT_REVIEW_SHA256).lower()
    review_raw = read_bytes(review_path, review_sha)
    review_value = json.loads(review_raw.decode("utf-8"))
    require(review_value.get("passed") is True, "Saved parent review did not pass")
    review = {
        "path": wire_path(review_path),
        "sha256": sha_bytes(review_raw),
        "bytes": len(review_raw),
        "schema": review_value.get("schema"),
        "component": review_value.get("component"),
        "passed": True,
        "outside_science_pins": True,
    }
    sympy = sympy_raw_pins(args.sympy_module_file or list(SYMPY_MODULE_DEFAULTS))
    kernel_path = resolve_path(args.kernel_path or REPO / KERNEL_RELATIVE, REPO)
    runner_path = resolve_path(args.runner_path or REPO / RUNNER_RELATIVE, REPO)
    require(kernel_path.resolve() == (REPO / KERNEL_RELATIVE).resolve(), "Kernel path escaped scope")
    require(runner_path.resolve() == (REPO / RUNNER_RELATIVE).resolve(), "Runner path escaped scope")
    require(kernel_path.exists() and runner_path.exists(), "Kernel and runner must exist before prepare")
    require(kernel_path.resolve() != Path(__file__).resolve(), "Metadata helper cannot be the kernel")
    if args.kernel_sha256:
        require(sha_path(kernel_path) == args.kernel_sha256.lower(), "Kernel SHA does not match")
    if args.runner_sha256:
        require(sha_path(runner_path) == args.runner_sha256.lower(), "Runner SHA does not match")
    input_raw, input_value = build_input_spec(
        parent,
        facts_path,
        facts,
        facts_raw,
        design_review,
        design_review_raw,
        review,
        sympy,
    )
    source = build_source(
        parent,
        facts_path,
        input_raw,
        facts,
        facts_raw,
        design_review,
        design_review_raw,
        review,
        sympy,
        kernel_path,
        runner_path,
    )
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
        "kernel": {"path": wire_path(kernel_path), "sha256": source["kernel_sha256"]},
        "runner": {"path": wire_path(runner_path), "sha256": source["runner_sha256"]},
        "metadata_only": True,
        "source_only": True,
        "no_runner_execution": True,
        "do_not_reenter": True,
    }
    require(len(encode(brief)) < 100_000, "Stage brief exceeds 100 KB")
    brief_raw = write_new(BRIEF, brief)
    pending_parent = {
        "schema": COMPONENT + "_pending_parent_authority_v1",
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "phase": "prepared_unentered",
        "metadata_only": True,
        "source_only": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
    }
    main_brief_raw, main_brief_value = update_main_brief(
        {
            "phase": "prepared_unentered",
            "component": COMPONENT,
            "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
            "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
            "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
            "source_union_count": source["source_union_count"],
            "parent_report": parent["report_meta"],
            "parent_result": parent["result_meta"],
            "pending_authoritative_parent": pending_parent,
            "no_math_imported": True,
            "no_math_called": True,
            "runner_executed": False,
            "metadata_only": True,
            "source_only": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
        }
    )
    receipt = {
        "schema": COMPONENT + "_prepare_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "input_spec_sha256": sha_bytes(input_raw),
        "diagnostic_checkpoint_sha256": sha_bytes(diagnostic_raw),
        "brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief_sha256": sha_bytes(main_brief_raw),
        "main_resume_brief": main_brief_value,
        "pending_authoritative_parent": pending_parent,
        "source_union_count": source["source_union_count"],
        "checkpoint_count": CURRENT_CHECKPOINT_COUNT,
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
        "source_raw_sha256": sha_bytes(source_raw),
        "kernel_path": source["kernel_path"],
        "kernel_sha256": source["kernel_sha256"],
        "derive_contract": "derive(source, diagnostic=None)",
        "derive_entry_calls": 0,
        "derive_calls": 0,
        "kernel_import_calls": 0,
        "metadata_only": True,
        "source_only": True,
        "runner_executed": False,
        "no_math_reexecuted": True,
    }
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    derive_entry_calls = 0
    derive_calls = 0
    kernel_import_calls = 0
    kernel_imported = False
    append_diagnostic(
        {
            "event": "manager_entry_started",
            "schema": DIAGNOSTIC_SCHEMA,
            "component": COMPONENT,
            "source_checkpoint_sha256": sha_bytes(source_raw),
            "entered_checkpoint_sha256": sha_bytes(entered_raw),
        }
    )
    try:
        module = import_kernel(kernel_path)
        kernel_import_calls = 1
        kernel_imported = True
        append_diagnostic(
            {
                "event": "manager_kernel_imported",
                "schema": DIAGNOSTIC_SCHEMA,
                "component": COMPONENT,
                "kernel_path": source["kernel_path"],
                "kernel_sha256": source["kernel_sha256"],
            }
        )

        def diagnostic(event: dict) -> None:
            # Bare kernel payloads are deliberately not wrapped or rewritten.
            append_diagnostic(event)

        derive_entry_calls = 1
        derive_calls = 1
        result = module.derive(source, diagnostic)
    except BaseException as error:  # persist the failed attempt before validation
        append_diagnostic(
            {
                "event": "manager_kernel_or_derive_failed",
                "schema": DIAGNOSTIC_SCHEMA,
                "component": COMPONENT,
                "error_type": type(error).__name__,
            }
        )
        result = failure_result(error, derive_entry_calls, derive_calls, kernel_import_calls)
    require(isinstance(result, dict), "Kernel result must be a JSON object")
    # Keep the kernel object exactly as a nested component result while the
    # durable outer envelope records manager-measured lifecycle facts.
    component_result = clone(result)
    durable_result = clone(result)
    durable_result["component_result"] = component_result
    durable_result["component_result_exact"] = True
    durable_result["execution_metadata"] = {
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_imported,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "new_math_entry": derive_calls,
        "old_math_reexecuted": False,
        "math_reexecuted": False,
        "runner_executed": False,
    }
    durable_result["kernel_import_calls"] = kernel_import_calls
    durable_result["kernel_imported"] = kernel_imported
    durable_result["derive_entry_calls"] = derive_entry_calls
    durable_result["derive_calls"] = derive_calls
    durable_result["new_math_entry"] = derive_calls
    durable_result["old_math_reexecuted"] = False
    durable_result["math_reexecuted"] = False
    durable_result["runner_executed"] = False
    result_raw = write_new(RESULT_CHECKPOINT, durable_result)
    diagnostic = diagnostic_summary()
    validation = validate_result(
        durable_result,
        source,
        derive_entry_calls,
        derive_calls,
        diagnostic["callback_event_count"],
        diagnostic["manager_kernel_imported"],
        diagnostic["callback_event_counts"],
        diagnostic["row_callback_event_count"],
    )
    receipt = {
        "schema": COMPONENT + "_run_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_path(ENTERED_CHECKPOINT),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "diagnostic_summary": diagnostic,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_import_calls": kernel_import_calls,
        "kernel_imported": kernel_imported,
        "validation": validation,
        "raw_result_saved_before_validation": True,
        "component_result_preserved": True,
        "rerun": False,
        "runner_executed": False,
        "runtime_model_metadata": result.get("runtime_model_metadata"),
    }
    receipt_raw = write_new(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw)}


def checkpoint_artifacts() -> list[dict]:
    return [artifact_entry(path, embed=True) for path in CHECKPOINTS]


def failed_attempt_artifacts(
    code2: list[dict],
    input_artifact: dict,
    current_checkpoints: list[dict],
) -> dict:
    """Return the exact seven artifacts owned by a failed attempt.

    A failed publication preserves the parent ledger and adds only the two
    scoped code files, the frozen input, and the four current checkpoints.
    This is intentionally disjoint from the source-facts and review pins.
    """
    additions = {entry["path"]: entry["sha256"] for entry in code2}
    additions[input_artifact["path"]] = input_artifact["sha256"]
    additions.update({entry["path"]: entry["sha256"] for entry in current_checkpoints})
    require(len(additions) == 7, "Failed attempt must add exactly seven artifacts")
    return additions


def alias_trace_pointer(alias_trace: dict) -> dict:
    return {
        "path": alias_trace["path"],
        "sha256": alias_trace["sha256"],
        "bytes": alias_trace["bytes"],
        "outside_science_pins": True,
        "static_only": True,
        "no_cross_side_parameter_identity_claim": True,
        "no_zero_or_error_cap_claim": True,
    }


def next_task_metadata(passed: bool, alias_trace: dict | None = None) -> dict:
    if not passed:
        return {
            "phase": "cached_W_sector_prefix_inspection",
            "action": "Inspect the saved W-sector result and diagnostic prefix only; do not reenter.",
            "source": "saved_W_sector_result_and_diagnostic_prefix",
            "reentry": "forbidden",
            "physical_layer_completed": False,
            "claims": "no zero, error-cap deletion, global closure, temporal recursion, or physical completion",
        }
    require(isinstance(alias_trace, dict), "Passed next task requires the alias trace pointer")
    pointer = alias_trace_pointer(alias_trace)
    return {
        "phase": "exact_source_parameter_pullback_for_rebuilt_W",
        "action": (
            "Define and bind source-native parameter interfaces for rebuilt W; distinguish its "
            "interval-origin delta enclosure from the actual native scalar/function."
        ),
        "source": "saved_W_sector_result_and_static_parameter_alias_trace",
        "parameter_alias_trace_pointer": pointer,
        "retain_unresolved": ["Delta R", "Delta S", "Delta delta", "pressure sectors", "history sectors"],
        "reentry": "forbidden",
        "physical_layer_completed": False,
        "claims": "no same-source parameter identity, zero, error-cap deletion, global closure, temporal recursion, or physical completion",
    }


def report_for_publish(
    source: dict,
    source_raw: bytes,
    result: dict,
    result_raw: bytes,
    validation: dict,
    code2: list[dict],
    input_artifact: dict,
    review_artifact: dict,
) -> dict:
    current = checkpoint_artifacts()
    checkpoint_hashes = {item["path"]: item["sha256"] for item in current}
    if validation["passed"]:
        failed_hashes = clone(source["parent_failed_hashes"])
    else:
        failed_hashes = merge_hashes(
            source["parent_failed_hashes"],
            failed_attempt_artifacts(code2, input_artifact, current),
        )
    component_result = result.get("component_result", result)
    helper = artifact_entry(Path(__file__).resolve(), embed=True, outside_science_pins=True)
    alias_trace = parameter_alias_trace_artifact()
    alias_pointer = alias_trace_pointer(alias_trace)
    next_task = next_task_metadata(validation["passed"], alias_trace if validation["passed"] else None)
    facts_path = resolve_path(source["sector_source_facts_path"], WORK)
    facts = artifact_entry(
        facts_path,
        SOURCE_FACTS_SHA256,
        embed=True,
        outside_science_pins=True,
    )
    facts_v1_path = resolve_path(source["sector_source_facts_v1_artifact"]["path"], WORK)
    facts_v1 = artifact_entry(
        facts_v1_path,
        SOURCE_FACTS_V1_SHA256,
        embed=True,
        outside_science_pins=True,
    )
    design_review = artifact_entry(
        DESIGN_REVIEW_DEFAULT,
        DESIGN_REVIEW_SHA256,
        embed=True,
        outside_science_pins=True,
    )
    sympy_artifacts = [
        {**entry, "raw_bytes_base64": entry["raw_bytes_base64"]}
        for entry in source["sympy_raw_pins"]["artifacts"]
    ]
    lossless = [
        *code2,
        input_artifact,
        facts,
        facts_v1,
        design_review,
        *current,
        helper,
        review_artifact,
        alias_trace,
        *sympy_artifacts,
    ]
    return {
        "schema": PUBLICATION_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if validation["passed"] else "failed",
        "passed": validation["passed"],
        "explicit_result_passed": result.get("passed", result.get("status") == "ok"),
        "result_status": result.get("status"),
        "metadata_only": True,
        "publication_metadata_only": True,
        "source_only": True,
        "conditional_source_functions": True,
        "result_metadata_only": result.get("metadata_only"),
        "result_source_only": result.get("source_only"),
        "component_result_metadata_only": component_result.get("metadata_only") if isinstance(component_result, dict) else None,
        "component_result_source_only": component_result.get("source_only") if isinstance(component_result, dict) else None,
        "physical_layer_completed": False,
        "pointwise_identity_proved": False,
        "cap_error_reduction_claimed": False,
        "accepted_progress_physical_baseline_preserved": True,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_parent_commit": PARENT_PARENT_COMMIT,
        "parent_report": source["parent_report"],
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
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
        "checkpoint_count": CURRENT_CHECKPOINT_COUNT,
        "checkpoint_paths": list(checkpoint_hashes),
        "checkpoint_hashes": checkpoint_hashes,
        "source_pin_count": source["source_union_count"],
        "source_union_count": source["source_union_count"],
        "source_pins": source["source_union_hashes"],
        "source_union_hashes": source["source_union_hashes"],
        "raw_sha_pins": source["raw_sha_pins"],
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
        "sector_source_facts_artifact": facts,
        "sector_source_facts_v1_artifact": facts_v1,
        "design_review_artifact": design_review,
        "kernel_contract": source["kernel_contract"],
        "sympy_raw_pins": source["sympy_raw_pins"],
        "input_artifact": input_artifact,
        "input_raw_sha256": input_artifact["sha256"],
        "code2": code2,
        "current_code_hashes": {entry["path"]: entry["sha256"] for entry in code2},
        "metadata_helper_artifact": helper,
        "external_draft": helper,
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
        "result_checkpoint_links": {
            "source": sha_bytes(source_raw),
            "entered": sha_path(ENTERED_CHECKPOINT),
            "result": sha_bytes(result_raw),
            "diagnostics": sha_path(DIAGNOSTIC_CHECKPOINT),
        },
        "actual_counts": validation["actual_counts"],
        "new_flags": validation["new_flags"],
        "protocol_violations": validation["protocol_violations"],
        "durable_result_validation": validation,
        "execution_metadata": clone(result.get("execution_metadata")),
        "runtime_model_metadata": result.get("runtime_model_metadata"),
        "kernel_import_calls": result.get("kernel_import_calls"),
        "kernel_imported": result.get("kernel_imported"),
        "derive_entry_calls": validation["derive_entry_calls"],
        "derive_calls": validation["derive_calls"],
        "new_math_entry": validation["derive_entry_calls"],
        "old_math_reexecuted": False,
        "math_reexecuted": False,
        "runner_executed": False,
        "no_runner_execution": True,
        "full_goal_complete": False,
        "goal_active": True,
        "do_not_reenter": True,
        "exact_nine_scope": True,
        "scoped_paths": [path.as_posix() for path in SCOPED_RELATIVES],
        "input_spec_schema": INPUT_SCHEMA,
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": list(NEW_FLAGS),
        "lossless_artifacts": lossless,
        "lossless_artifact_count": len(lossless),
        "independent_review": review_artifact,
        "parameter_alias_trace_artifact": alias_trace,
        "parameter_alias_trace_pointer": alias_pointer,
        "next_executable_task": next_task,
    }


def update_main_brief(patch: dict) -> tuple[bytes, dict]:
    require(MAIN_BRIEF.exists(), "Main resume brief is missing")
    old_raw, old = read_json(MAIN_BRIEF)
    require(len(old_raw) < 100_000, "Main resume brief exceeds 100 KB")
    preimage = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(old_raw)[:12]}.json"
    if not preimage.exists():
        write_new(preimage, old_raw)
    else:
        require(sha_path(preimage) == sha_bytes(old_raw), "Existing BRIEF preimage changed")
    result = brief_safe(clone(old))
    in_progress = result.get("in_progress")
    if not isinstance(in_progress, dict):
        in_progress = {}
    in_progress[COMPONENT] = brief_safe(portable(patch))
    result["in_progress"] = in_progress
    result["latest_component"] = COMPONENT
    result["latest_component_pointer"] = brief_safe(portable(patch))
    result["current_exact_source_result"] = brief_safe(portable(patch))
    canonical = [
        wire_path(REPO / "AGENTS.md"),
        wire_path(REPO / MILESTONE_RELATIVE),
        wire_path(REPO / REPORT_RELATIVE),
    ]
    old_read_first = result.get("read_first")
    if not isinstance(old_read_first, list):
        old_read_first = []
    published = (REPO / REPORT_RELATIVE).exists() and (REPO / MILESTONE_RELATIVE).exists()
    if published:
        result["read_first"] = canonical
        result["read_first_alias_repairs"] = {
            "AGENTS.md": canonical[0],
            "newdoc": canonical[1],
            "report": canonical[2],
        }
    else:
        result["read_first"] = clone(old_read_first)
    for alias in (
        "latest_review",
        "latest_postrun_saved_review",
        "latest_qa",
        "latest_final_qa",
        "latest_seal",
        "latest_final_seal",
        "current_document_sha256",
        "current_doc_sha256",
        "next_executable_task",
        "next_executable_task_source",
        "next_remaining",
        "next_remaining_scope_plan",
    ):
        if alias in patch:
            result[alias] = brief_safe(portable(patch[alias]))
    authority = patch.get("latest_authoritative_state")
    if authority is not None:
        require(isinstance(authority, dict), "Authoritative state is not an object")
        safe_authority = brief_safe(portable(authority))
        result["latest_authoritative_state"] = safe_authority
        result["verified_local_and_remote_commit"] = safe_authority.get(
            "verified_local_and_remote_commit"
        )
        if "next_remaining" in safe_authority:
            result["next_remaining"] = safe_authority["next_remaining"]
    if published:
        milestone_sha256 = sha_path(REPO / MILESTONE_RELATIVE)
        result["latest_milestone"] = {
            "path": canonical[1],
            "report": canonical[2],
            "phase": patch.get("phase"),
            "report_sha256": patch.get("report_sha256"),
            "milestone_sha256": milestone_sha256,
        }
    raw = encode(result)
    if len(raw) > 100_000:
        compact_keys = {
            "continuation_detail_archive",
            "archived_detail",
            "archived_prior_component_details",
            "failed_source_attempt_evidence",
        }
        preserve_keys = (
            "latest_component",
            "latest_component_pointer",
            "current_exact_source_result",
            "execution_rules",
            "legacy_dirty_paths",
            "legacy_dirty_path_hashes",
            "protected_file_hashes",
            "next_executable_task",
            "next_remaining",
            "next_remaining_scope_plan",
        )
        preserved_values = {key: result.get(key) for key in preserve_keys if key in result}
        previous_in_progress = result.get("in_progress")
        current_in_progress = None
        previous_selectors = {}
        if isinstance(previous_in_progress, dict):
            current_in_progress = previous_in_progress.get(COMPONENT)
            previous_selectors = {
                key: value
                for key, value in previous_in_progress.items()
                if key != COMPONENT
            }
        compact = {key: value for key, value in result.items() if key not in compact_keys}
        if isinstance(previous_in_progress, dict):
            selector_pointer = {
                "path": wire_path(preimage),
                "sha256": sha_bytes(old_raw),
                "bytes": len(old_raw),
                "preserved_exact_bytes": True,
            }
            compact_in_progress = {}
            if current_in_progress is not None:
                compact_in_progress[COMPONENT] = current_in_progress
            for selector in sorted(previous_selectors):
                compact_in_progress[selector] = {
                    "archived_to_main_resume_brief_preimage": selector_pointer,
                    "selector": selector,
                }
            compact["in_progress"] = compact_in_progress
            compact["archived_in_progress_selectors"] = {
                "preimage": selector_pointer,
                "selectors": sorted(previous_selectors),
            }
        compact["archived_main_resume_brief_preimage"] = {
            "path": wire_path(preimage),
            "sha256": sha_bytes(old_raw),
            "bytes": len(old_raw),
            "preserved_exact_bytes": True,
        }
        compact["main_resume_brief_archival_compaction"] = True
        result = compact
        raw = encode(result)
        for key, value in preserved_values.items():
            require(result.get(key) == value, f"BRIEF compaction dropped protected metadata: {key}")
    require(len(raw) < 100_000, "Updated main resume brief exceeds 100 KB")
    write_replace(MAIN_BRIEF, raw)
    return raw, result


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
    require(valid_git_revision(head), "Local Git head is not a 40-hex revision")
    require(git_output("rev-parse", "HEAD^").lower() == PARENT_COMMIT, "Published parent changed")
    branch = git_output("branch", "--show-current")
    require(branch, "Published branch is unnamed")
    remote = git_output("ls-remote", "origin", "refs/heads/" + branch).split()
    require(remote and valid_git_revision(remote[0]), "Remote branch head is unavailable")
    require(remote[0].lower() == head.lower(), "Published head is not synchronized")
    changed = set(git_output("diff", "--name-only", PARENT_COMMIT + "..HEAD").splitlines())
    require(changed == {path.as_posix() for path in SCOPED_RELATIVES}, "Committed scope is not exactly nine paths")
    validate_scope_attributes()
    for relative in SCOPED_RELATIVES:
        local_path = REPO / relative
        require(local_path.exists(), "Published scoped path is missing: " + relative.as_posix())
        blob_size = int(git_output("cat-file", "-s", "HEAD:" + relative.as_posix()))
        require(blob_size == local_path.stat().st_size, "Committed blob byte count changed")
        blob = git_output("rev-parse", "HEAD:" + relative.as_posix()).lower()
        require(valid_git_revision(blob), "Committed blob is not a Git object")
        require(blob == git_output("hash-object", wire_path(local_path)).lower(), "Published bytes differ")
    require(report.get("parent_commit") == PARENT_COMMIT, "Report parent commit changed")
    return head


def publish(args: argparse.Namespace) -> dict:
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT)
    _, run_receipt = read_json(RUN_RECEIPT)
    validation = validate_result(
        result,
        source,
        run_receipt.get("derive_entry_calls"),
        run_receipt.get("derive_calls"),
        run_receipt.get("diagnostic_summary", {}).get("callback_event_count"),
        run_receipt.get("diagnostic_summary", {}).get("manager_kernel_imported"),
        run_receipt.get("diagnostic_summary", {}).get("callback_event_counts"),
        run_receipt.get("diagnostic_summary", {}).get("row_callback_event_count"),
    )
    kernel = resolve_path(source["kernel_path"], REPO)
    runner = resolve_path(source["runner_path"], REPO)
    code2 = [
        artifact_entry(kernel, source["kernel_sha256"], embed=True),
        artifact_entry(runner, source["runner_sha256"], embed=True),
    ]
    input_artifact = artifact_entry(
        resolve_path(source["input_spec_metadata"]["path"], WORK),
        source["input_spec_metadata"]["sha256"],
        embed=True,
    )
    review_artifact = artifact_entry(
        resolve_path(source["review"]["path"], WORK),
        source["review"]["sha256"],
        embed=True,
        outside_science_pins=True,
    )
    report = report_for_publish(
        source, source_raw, result, result_raw, validation, code2, input_artifact, review_artifact
    )
    report_raw = write_replace(REPO / REPORT_RELATIVE, report)
    milestone = (
        "# Current original O3 transition W sector correlation (2026-10-27)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Status: `{report['attempt_status']}`; actual counts: `{report['actual_counts']}`.\n\n"
        "This is a metadata-only publication envelope carrying a source-conditional "
        "full-W symbolic sector leaf. It records "
        "eight signed sector families, 32 numerator jet rows, four cross-numerator "
        "identity orders, and four unevaluated views of D_Z^j[(sum of eight cross "
        "numerator sectors)/(L_n L_c)]. It does not claim a physical point identity, "
        "P2/W zero, cap reduction, numeric containment, global closure, temporal "
        "recursion, phase equality, or physical completion.\n\n"
        "The saved result is immutable and must not be reentered.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(REPO / MILESTONE_RELATIVE, milestone)
    prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; do not reenter.\n"
        f"Next phase: `{report['next_executable_task']['phase']}`.\n"
        f"Next action: {report['next_executable_task']['action']}\n"
        f"Next action source: `{report['next_executable_task']['source']}`.\n"
        f"Static parameter-alias trace (outside science pins): `{report['parameter_alias_trace_pointer']['path']}` "
        f"SHA-256 `{report['parameter_alias_trace_pointer']['sha256']}`; no cross-side identity, zero, or error-cap claim.\n"
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
    brief_raw, brief = update_main_brief(
        {
            "phase": "published_passed" if report["passed"] else "published_failed",
            "report": wire_path(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "milestone": wire_path(REPO / MILESTONE_RELATIVE),
            "milestone_sha256": sha_bytes(milestone_raw),
            "attempt_status": report["attempt_status"],
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "latest_review": review_artifact,
            "latest_postrun_saved_review": review_artifact,
            "latest_qa": validation,
            "latest_final_qa": validation,
            "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
            "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
            "current_document_sha256": sha_bytes(milestone_raw),
            "current_doc_sha256": sha_bytes(milestone_raw),
            "next_executable_task": report.get("next_executable_task"),
            "next_executable_task_source": report.get("next_executable_task", {}).get("source"),
            "next_remaining": report.get("next_executable_task"),
            "next_remaining_scope_plan": report.get("next_executable_task"),
        }
    )
    receipt = {
        "schema": COMPONENT + "_publication_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "handoff_preimages": preimages,
        "exact_nine_scope": True,
        "runtime_model_metadata": None,
        "main_resume_brief": {"sha256": sha_bytes(brief_raw), "value": brief},
    }
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
    validation = validate_result(
        result,
        source,
        stored.get("derive_entry_calls"),
        stored.get("derive_calls"),
        stored.get("diagnostic_callback_count"),
        stored.get("diagnostic_manager_kernel_imported"),
        stored.get("diagnostic_event_counts"),
        stored.get("diagnostic_row_callback_count"),
    )
    require(validation == stored, "Durable validation changed")
    head = verify_published_commit(report)
    milestone_sha256 = sha_path(REPO / MILESTONE_RELATIVE)
    brief_raw, brief = update_main_brief(
        {
            "phase": "resumed_published",
            "report": wire_path(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "head": head,
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "read_first": [
                wire_path(REPO / "AGENTS.md"),
                wire_path(REPO / MILESTONE_RELATIVE),
                wire_path(REPO / REPORT_RELATIVE),
            ],
            "latest_review": report.get("independent_review"),
            "latest_postrun_saved_review": report.get("independent_review"),
            "latest_qa": validation,
            "latest_final_qa": validation,
            "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
            "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "pending"},
            "current_document_sha256": milestone_sha256,
            "current_doc_sha256": milestone_sha256,
            "next_executable_task": report.get("next_executable_task"),
            "next_executable_task_source": report.get("next_executable_task", {}).get("source"),
            "next_remaining": report.get("next_executable_task"),
            "next_remaining_scope_plan": report.get("next_executable_task"),
        }
    )
    receipt = {
        "schema": COMPONENT + "_resume_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "head": head,
        "main_resume_brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief": brief,
        "do_not_reenter": True,
    }
    receipt_raw = write_replace(RESUME_RECEIPT, receipt)
    return {"resume": receipt, "resume_sha256": sha_bytes(receipt_raw)}


def seal(args: argparse.Namespace) -> dict:
    review_path = resolve_path(args.review or PARENT_REVIEW_DEFAULT, WORK)
    review_sha = (args.review_sha256 or PARENT_REVIEW_SHA256).lower()
    review_raw = read_bytes(review_path, review_sha)
    review = json.loads(review_raw.decode("utf-8"))
    require(review.get("passed") is True, "Independent review did not pass")
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    stored = report.get("durable_result_validation")
    require(isinstance(stored, dict), "Durable validation is missing before seal")
    validation = validate_result(
        result,
        source,
        stored.get("derive_entry_calls"),
        stored.get("derive_calls"),
        stored.get("diagnostic_callback_count"),
        stored.get("diagnostic_manager_kernel_imported"),
        stored.get("diagnostic_event_counts"),
        stored.get("diagnostic_row_callback_count"),
    )
    require(validation == stored, "Durable validation changed before seal")
    require(report.get("full_result") == result, "Published result changed before seal")
    head = verify_published_commit(report)
    milestone_sha256 = sha_path(REPO / MILESTONE_RELATIVE)
    authority = {
        "schema": COMPONENT + "_authoritative_state_v1",
        "component": COMPONENT,
        "commit": head,
        "head": head,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "passed": report.get("passed") is True,
        "attempt_status": report.get("attempt_status"),
        "actual_counts": clone(report.get("actual_counts")),
        "new_flags": clone(report.get("new_flags")),
        "full_goal_complete": False,
        "physical_layer_completed": report.get("physical_layer_completed") is True,
        "metadata_only": report.get("metadata_only") is True,
        "source_only": report.get("source_only") is True,
        "verified_local_and_remote_commit": True,
        "do_not_reenter": True,
        "review": {"path": wire_path(review_path), "sha256": sha_bytes(review_raw), "bytes": len(review_raw)},
        "next_remaining": report.get("next_executable_task"),
        "read_first": [
            wire_path(REPO / "AGENTS.md"),
            wire_path(REPO / MILESTONE_RELATIVE),
            wire_path(REPO / REPORT_RELATIVE),
        ],
    }
    brief_raw, brief = update_main_brief({
        "phase": "sealed",
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "seal": wire_path(SEAL_RECEIPT),
        "head": head,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "latest_authoritative_state": authority,
        "latest_review": authority["review"],
        "latest_postrun_saved_review": authority["review"],
        "latest_qa": validation,
        "latest_final_qa": validation,
        "latest_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "sealed"},
        "latest_final_seal": {"path": wire_path(SEAL_RECEIPT), "phase": "sealed"},
        "current_document_sha256": milestone_sha256,
        "current_doc_sha256": milestone_sha256,
        "next_executable_task": report.get("next_executable_task"),
        "next_executable_task_source": report.get("next_executable_task", {}).get("source"),
        "next_remaining": report.get("next_executable_task"),
        "next_remaining_scope_plan": report.get("next_executable_task"),
    })
    receipt = {
        "schema": SEAL_SCHEMA,
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "review": authority["review"],
        "head": head,
        "main_resume_brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief": brief,
        "do_not_reenter": True,
    }
    receipt_raw = write_replace(SEAL_RECEIPT, receipt)
    return {"seal": receipt, "seal_sha256": sha_bytes(receipt_raw)}


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("action", choices=("prepare", "run", "publish", "resume", "seal"))
    value.add_argument("--parent-report")
    value.add_argument("--parent-report-sha256")
    value.add_argument("--source-facts")
    value.add_argument("--source-facts-sha256")
    value.add_argument("--review")
    value.add_argument("--review-sha256")
    value.add_argument("--kernel-path")
    value.add_argument("--runner-path")
    value.add_argument("--kernel-sha256")
    value.add_argument("--runner-sha256")
    value.add_argument("--sympy-module-file", action="append", default=[])
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
