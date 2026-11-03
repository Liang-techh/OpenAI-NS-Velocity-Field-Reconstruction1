"""Run the scoped graph evaluator boundary exactly once.

This runner has a deliberately narrow boundary.  It authenticates the raw
leaf, fixture, runner, graph archive, and prior program receipt before any
gzip decode.  It writes the exclusive entry record before loading the new
leaf or fixture, runs the synthetic fixture exactly once, and only then
decodes the two pinned archives.

The post-fixture work is a structural catalog only.  It reads the seven
``definite_integral`` ``original_record`` values from
``receipt['result']['program']['instructions']`` and scans the saved
``exact_graph_nodes`` list for same-name ``bound_variable`` rows.  It does
not follow graph edges, install binder mappings, call ``evaluate`` on the
real graph, import an old producer, or compute a source value.  A variable
name with zero or multiple candidates remains unresolved/ambiguous rather
than being guessed.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "n1_scoped_graph_evaluator_v1"
PHASE = COMPONENT + "_pending"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_scoped_graph_evaluator_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_scoped_graph_evaluator_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_scoped_graph_evaluator_v1_once.py"
GRAPH_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
PRIOR_REL = "experiments/root_st073/n1_raw_g5_dag_program_v1_receipt.json.gz"

FROZEN_INPUTS = (LEAF_REL, FIXTURE_REL, RUNNER_REL, GRAPH_REL, PRIOR_REL)
FROZEN_SET = frozenset(FROZEN_INPUTS)

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
PRIOR_SHA256 = "fa7108686ed6d6331bb8f09a12d86dba5ccd821d567aa45fbb1efc301f6f657f"

LEAF_SCHEMA = "n1_scoped_graph_evaluator_v1"
FIXTURE_SCHEMA = "n1_scoped_graph_evaluator_fixture_v1"
PRIOR_SCHEMA = "n1_raw_g5_dag_program_v1"
PRIOR_COMPONENT = "n1_raw_g5_dag_program_v1"

SHA256_HEX = frozenset("0123456789abcdefABCDEF")
_MISSING = object()


class RunnerInputError(ValueError):
    """Raised for a fail-closed input, provenance, or scope violation."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    """Serialize ordinary receipt data without invoking graph operations."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "opaque_bytes", "length": len(value), "sha256": _sha_bytes(value)}

    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "opaque_cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        return {"type": "opaque_value", "class": type(value).__name__, "repr": repr(value)}
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    descriptor = os.open(str(path), flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    finally:
        if descriptor != -1:
            os.close(descriptor)
    return _sha_bytes(raw)


def _decode_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _sha_text(value: Any, label: str) -> str:
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in SHA256_HEX for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " path must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " path is missing or escapes repository: " + key,
    )
    return key, path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list), "source.files list is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(files):
        _require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_file(row.get("path"), f"source.files[{index}]")
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        expected = _sha_text(row.get("sha256"), key)
        actual = _sha_file(path)
        _require(actual == expected, "raw SHA mismatch: " + key)
        if key in (GRAPH_REL, PRIOR_REL):
            _require(
                row.get("raw_sha256_before_decode") is True,
                "gzip input is not marked raw-pinned before decode: " + key,
            )
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET, "manifest must contain exactly the five frozen inputs")
    _require(hashes[GRAPH_REL] == GRAPH_SHA256, "real graph raw SHA is not the accepted pin")
    _require(hashes[PRIOR_REL] == PRIOR_SHA256, "prior receipt raw SHA is not the accepted pin")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _validate_contract(source: Mapping[str, Any]) -> None:
    contract = source.get("input_contract")
    if contract is None:
        return
    _require(isinstance(contract, Mapping), "input_contract must be an object")

    if "fixture_api" in contract:
        _require(contract.get("fixture_api") == "run_fixture(leaf)", "fixture API contract mismatch")
    if "graph_path" in contract:
        _require(contract.get("graph_path") == GRAPH_REL, "graph path contract mismatch")
    if "real_graph_path" in contract:
        _require(contract.get("real_graph_path") == GRAPH_REL, "real graph path contract mismatch")
    if "prior_receipt_path" in contract:
        _require(contract.get("prior_receipt_path") == PRIOR_REL, "prior receipt path contract mismatch")
    if "prior_program_path" in contract:
        _require(
            contract.get("prior_program_path") == ["result", "program", "instructions"],
            "prior program path contract mismatch",
        )
    for key in (
        "authenticate_raw_gzip_before_decode",
        "decode_after_fixture",
        "no_graph_import",
        "no_graph_traversal",
        "no_graph_edge_traversal",
        "no_actual_graph_evaluation",
        "no_source_evaluation",
    ):
        if key in contract:
            _require(contract.get(key) is True, key + " contract is disabled")
    if "integral_record_count" in contract:
        _require(contract.get("integral_record_count") == 7, "integral count contract mismatch")


def _brief_record(brief: Mapping[str, Any]) -> Mapping[str, Any]:
    if brief.get("component") == COMPONENT:
        return brief
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        return nested
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping) and continuation.get("component") == COMPONENT:
        return continuation
    return brief


def _checkpoint_candidate(brief: Mapping[str, Any], record: Mapping[str, Any]) -> Mapping[str, Any] | None:
    candidates: list[Any] = [record.get("source_checkpoint"), brief.get("source_checkpoint")]
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        candidates.append(nested.get("source_checkpoint"))
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping):
        candidates.append(continuation.get("source_checkpoint"))
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            return candidate
    return None


def _checkpoint_path(raw_path: Any) -> Path:
    _require(isinstance(raw_path, str) and raw_path, "pending source checkpoint path is missing")
    path_value = Path(raw_path.replace("\\", "/"))
    if path_value.is_absolute():
        return path_value.resolve()
    _require(".." not in path_value.parts, "pending source checkpoint escapes repository")
    return (ROOT / path_value).resolve()


def _validate_pending_brief(source_sha: str) -> tuple[str, str]:
    _require(BRIEF.is_file(), "pending resume brief is missing")
    brief_raw = BRIEF.read_bytes()
    brief_sha = _sha_bytes(brief_raw)
    brief = _decode_object(brief_raw, "pending resume brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending brief one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending brief re-entry guard missing")
    if "goal_active" in record:
        _require(record.get("goal_active") is True, "pending brief goal is not active")
    if "full_goal_complete" in record:
        _require(record.get("full_goal_complete") is False, "pending brief marks the full goal complete")

    checkpoint = _checkpoint_candidate(brief, record)
    _require(checkpoint is not None, "pending brief source_checkpoint is missing")
    checkpoint_path = _checkpoint_path(checkpoint.get("path"))
    _require(checkpoint_path == SOURCE.resolve(), "pending brief points to another source checkpoint")
    _require(
        _sha_text(checkpoint.get("sha256"), "pending source checkpoint") == source_sha,
        "pending brief source checkpoint SHA mismatch",
    )
    return brief_sha, SOURCE.relative_to(ROOT).as_posix()


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered or has a result")
    _require(SOURCE.is_file(), "source manifest is missing")
    source_raw = SOURCE.read_bytes()
    source_sha = _sha_bytes(source_raw)
    source = _decode_object(source_raw, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("phase") == PHASE, "source manifest phase mismatch")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")
    _require(source.get("do_not_reenter") is True, "source manifest re-entry guard missing")
    _validate_contract(source)
    hashes = _manifest_hashes(source)
    brief_sha, checkpoint_key = _validate_pending_brief(source_sha)
    return source, hashes, source_sha, brief_sha, checkpoint_key


def _read_frozen_raw(hashes: Mapping[str, str]) -> dict[str, bytes]:
    """Read and re-authenticate bytes without decoding either gzip input."""

    raw: dict[str, bytes] = {}
    for relative in FROZEN_INPUTS:
        path = ROOT / relative
        data = path.read_bytes()
        _require(_sha_bytes(data) == hashes[relative], "frozen input changed after preflight: " + relative)
        raw[relative] = data
    return raw


def _forbidden_project_modules() -> list[str]:
    """Reject historical package/source modules already present in the process."""

    found: list[str] = []
    package_root = (ROOT / "src" / "openai_ns_reconstruction").resolve()
    experiments_root = (ROOT / "experiments").resolve()
    for name, module in list(sys.modules.items()):
        if name in _LOADED_NAMES:
            continue
        if name.startswith("openai_ns_reconstruction"):
            found.append(name)
            continue
        module_file = getattr(module, "__file__", None)
        if not isinstance(module_file, str):
            continue
        try:
            resolved = Path(module_file).resolve()
        except OSError:
            continue
        if resolved.is_relative_to(package_root) or resolved.is_relative_to(experiments_root):
            found.append(name)
    return sorted(set(found))


_LOADED_NAMES = (
    "_n1_scoped_graph_evaluator_v1_runner_leaf",
    "_n1_scoped_graph_evaluator_v1_runner_fixture",
)


def _cleanup_new_modules() -> None:
    for name in _LOADED_NAMES:
        sys.modules.pop(name, None)


def _load_module(relative: str, module_name: str) -> ModuleType:
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(module_name, path)
    _require(spec is not None and spec.loader is not None, "new module spec is unavailable: " + relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _load_new_modules() -> tuple[ModuleType, ModuleType]:
    _require(not _forbidden_project_modules(), "historical package or source module already loaded")
    try:
        leaf = _load_module(LEAF_REL, _LOADED_NAMES[0])
        fixture = _load_module(FIXTURE_REL, _LOADED_NAMES[1])
        _require(not _forbidden_project_modules(), "historical module loaded by new files")
        return leaf, fixture
    except BaseException:
        _cleanup_new_modules()
        raise


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not a mapping")
    result = dict(value)
    _require(result.get("schema") == FIXTURE_SCHEMA, "fixture schema mismatch")
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    _require(result.get("synthetic_only") is True, "fixture is not synthetic-only")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks), "fixture checks must be a non-empty mapping")
    _require(all(item is True for item in checks.values()), "synthetic fixture contains a failed check")
    for key in (
        "used_old_results",
        "used_old_module",
        "historical_module_imported",
        "historical_modules_imported",
        "old_receipt_decoded",
        "old_receipt_decompressed",
        "graph_decoded",
        "graph_decompressed",
        "source_graph_imported",
        "source_graph_executed",
        "actual_graph_evaluated",
    ):
        if key in result:
            _require(result[key] in (False, 0, None), "fixture reported forbidden activity: " + key)
    return _json_safe(result)


def _decode_gzip_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
        value = json.loads(decoded)
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid authenticated gzip JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _prior_program(receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(receipt.get("component") == PRIOR_COMPONENT, "prior receipt component mismatch")
    result = receipt.get("result")
    _require(isinstance(result, Mapping), "prior receipt result is missing")
    program = result.get("program")
    _require(isinstance(program, Mapping), "prior receipt result.program is missing")
    _require(program.get("schema") == PRIOR_SCHEMA, "prior program schema mismatch")
    instructions = program.get("instructions")
    _require(isinstance(instructions, list), "prior program instructions are missing")
    return program


def _integral_records(program: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Read only integral original records from the saved instruction list."""

    instructions = program.get("instructions")
    _require(isinstance(instructions, list), "prior program instructions are missing")
    records: list[dict[str, Any]] = []
    seen_ids: set[int] = set()
    for instruction_index, instruction in enumerate(instructions):
        _require(isinstance(instruction, Mapping), f"program instruction {instruction_index} is not an object")
        original = instruction.get("original_record", _MISSING)
        if not isinstance(original, Mapping):
            continue
        operation = original.get("operation", _MISSING)
        if operation != "definite_integral":
            continue
        node_id = instruction.get("node_id", _MISSING)
        _require(type(node_id) is int and node_id >= 0, "integral instruction node_id is malformed")
        _require(node_id not in seen_ids, "duplicate integral instruction node_id")
        seen_ids.add(node_id)
        variable = original.get("variable", _MISSING)
        _require(isinstance(variable, str) and variable, f"integral {node_id} variable name is missing")
        records.append(
            {
                "instruction_index": instruction_index,
                "integral_node_id": node_id,
                "variable_name": variable,
                "instruction_status": instruction.get("status"),
                "instruction_op": instruction.get("op"),
                "original_record": dict(original),
            }
        )
    return records


def _exact_graph_nodes(graph_archive: Mapping[str, Any]) -> list[Any]:
    nodes = graph_archive.get("exact_graph_nodes", _MISSING)
    _require(isinstance(nodes, list) and bool(nodes), "real graph exact_graph_nodes list is missing")
    return nodes


def _scan_bound_variable_candidates(
    graph_archive: Mapping[str, Any], variable_names: list[str]
) -> tuple[dict[str, list[int]], int]:
    """Scan row metadata only; never inspect an edge/reference field."""

    nodes = _exact_graph_nodes(graph_archive)
    wanted = set(variable_names)
    candidates = {name: [] for name in variable_names}
    for node_id, row in enumerate(nodes):
        if not isinstance(row, Mapping):
            continue
        if row.get("operation") != "bound_variable":
            continue
        name = row.get("name")
        if isinstance(name, str) and name in wanted:
            candidates[name].append(node_id)
    return candidates, len(nodes)


def _build_integral_catalog(
    graph_archive: Mapping[str, Any],
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    variable_names = list(dict.fromkeys(record["variable_name"] for record in records))
    candidates, node_count = _scan_bound_variable_candidates(graph_archive, variable_names)
    catalog_records: list[dict[str, Any]] = []
    unique_count = 0
    ambiguous_count = 0
    unresolved_count = 0
    for record in records:
        name = record["variable_name"]
        node_ids = list(candidates[name])
        if len(node_ids) == 1:
            status = "unique_candidate"
            unique_id: int | None = node_ids[0]
            unique_count += 1
        elif not node_ids:
            status = "unresolved_no_candidate"
            unique_id = None
            unresolved_count += 1
        else:
            status = "ambiguous_multiple_candidates"
            unique_id = None
            ambiguous_count += 1
        catalog_records.append(
            {
                **record,
                "candidate_node_ids": node_ids,
                "unique_node_id": unique_id,
                "binding_status": status,
                "binder_mapping_installed": False,
            }
        )

    by_variable = {
        name: {
            "candidate_node_ids": list(candidates[name]),
            "unique_node_id": candidates[name][0] if len(candidates[name]) == 1 else None,
            "binding_status": (
                "unique_candidate"
                if len(candidates[name]) == 1
                else "unresolved_no_candidate"
                if not candidates[name]
                else "ambiguous_multiple_candidates"
            ),
            "binder_mapping_installed": False,
        }
        for name in variable_names
    }
    return {
        "schema": COMPONENT + "_integral_catalog_v1",
        "integral_record_count": len(catalog_records),
        "expected_integral_record_count": 7,
        "records": catalog_records,
        "by_variable_name": by_variable,
        "variable_names": variable_names,
        "graph_node_count_scanned": node_count,
        "scan_operation": "bound_variable name metadata only",
        "graph_edges_traversed": False,
        "actual_graph_evaluate_called": False,
        "actual_source_values_computed": False,
        "binder_mappings_installed": False,
        "unique_candidate_count": unique_count,
        "ambiguous_candidate_count": ambiguous_count,
        "unresolved_candidate_count": unresolved_count,
        "unknown_duplicates_not_bound": True,
    }


def _consume(
    hashes: Mapping[str, str],
    progress: dict[str, Any],
) -> dict[str, Any]:
    # Raw bytes are authenticated before either gzip is decoded.  This also
    # detects a post-entry replacement before the new modules are loaded.
    raw = _read_frozen_raw(hashes)
    leaf, fixture = _load_new_modules()
    progress["new_modules_imported"] = True

    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "fixture.run_fixture is missing or not callable")
    progress["fixture_calls"] = 1
    try:
        raw_fixture = run_fixture(leaf)
    except BaseException:
        progress["fixture_traceback"] = traceback.format_exc()
        raise
    progress["fixture_raw"] = _json_safe(raw_fixture)
    progress["fixture"] = progress["fixture_raw"]
    fixture_result = _validate_fixture(raw_fixture)
    progress["fixture"] = fixture_result
    progress["synthetic_fixture_values_computed"] = True
    _require(not _forbidden_project_modules(), "historical module loaded during fixture")

    # Only after the fixture passes may either pinned gzip be decoded.
    raw = _read_frozen_raw(hashes)
    graph_archive = _decode_gzip_json(raw[GRAPH_REL], GRAPH_REL)
    progress["graph_decoded"] = True
    receipt = _decode_gzip_json(raw[PRIOR_REL], PRIOR_REL)
    progress["prior_receipt_decoded"] = True
    program = _prior_program(receipt)

    records = _integral_records(program)
    progress["integral_records_partial"] = _json_safe(records)
    _require(len(records) == 7, f"expected exactly seven integral original_records, found {len(records)}")

    catalog = _build_integral_catalog(graph_archive, records)
    progress["integral_catalog"] = _json_safe(catalog)

    return {
        "status": "returned",
        "passed": True,
        "fixture": fixture_result,
        "integral_catalog": catalog,
        "catalog": catalog,
        "prior_receipt_component": receipt.get("component"),
        "prior_program_schema": program.get("schema"),
        "prior_program_path": ["result", "program", "instructions"],
        "real_graph_path": GRAPH_REL,
        "real_graph_raw_sha256": hashes[GRAPH_REL],
        "prior_receipt_path": PRIOR_REL,
        "prior_receipt_raw_sha256": hashes[PRIOR_REL],
        "fixture_calls": 1,
        "actual_integral_records_extracted": 7,
        "actual_graph_evaluate_called": False,
        "actual_source_value_evaluated": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "graph_traversed": False,
        "graph_edges_traversed": False,
        "old_producer_imported": False,
        "graph_import_used": False,
        "numeric_proof": False,
        "certified": False,
        "synthetic_fixture_values_computed": True,
        "binder_mappings_installed": False,
        "unknown_duplicates_not_bound": True,
        "traceback": None,
    }


def _base_result(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    checkpoint_key: str,
    entry_sha: str,
    progress: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": COMPONENT + "_runner_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "source_checkpoint_path": SOURCE.relative_to(ROOT).as_posix(),
        "pending_brief_path": BRIEF.relative_to(ROOT.parent.parent).as_posix(),
        "pending_brief_sha256": brief_sha,
        "pending_brief_checkpoint_path": checkpoint_key,
        "source_checkpoint_sha256": source_sha,
        "input_hashes": dict(hashes),
        "one_entry_only": True,
        "do_not_reenter": True,
        "brief_modified": False,
        "raw_gzip_authenticated_before_decode": True,
        "decode_after_fixture_only": True,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "graph_traversed": False,
        "graph_edges_traversed": False,
        "old_producer_imported": False,
        "graph_import_used": False,
        "actual_graph_evaluate_called": False,
        "actual_source_value_evaluated": False,
        "synthetic_fixture_values_computed": False,
        "numeric_proof": False,
        "certified": False,
        "fixture_calls": int(progress.get("fixture_calls", 0)),
        "integral_record_count": len(progress.get("integral_records_partial", [])),
        "progress": _json_safe(progress),
        "source_manifest": _json_safe(source),
        "entry_sha256": entry_sha,
    }


def main() -> int:
    try:
        source, hashes, source_sha, brief_sha, checkpoint_key = _preflight()
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_runner_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed_preflight",
            "passed": False,
            "fixture": None,
            "integral_catalog": None,
            "traceback": traceback.format_exc(),
            "error": type(exc).__name__ + ": " + str(exc),
            "entry_written": False,
            "result_written": False,
        }
        print(json.dumps(output, ensure_ascii=False, sort_keys=True))
        return 1

    progress: dict[str, Any] = {
        "fixture_calls": 0,
        "new_modules_imported": False,
        "graph_decoded": False,
        "prior_receipt_decoded": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "graph_traversed": False,
        "graph_edges_traversed": False,
        "actual_graph_evaluate_called": False,
        "actual_source_value_evaluated": False,
        "integral_records_partial": [],
        "integral_catalog": None,
    }

    try:
        entry_sha = _write_exclusive(
            ENTRY,
            {
                "schema": COMPONENT + "_entered_v1",
                "component": COMPONENT,
                "phase": PHASE,
                "entered": True,
                "one_entry_only": True,
                "do_not_reenter": True,
                "source_checkpoint_sha256": source_sha,
                "pending_brief_sha256": brief_sha,
                "pending_brief_checkpoint_path": checkpoint_key,
                "input_hashes": dict(hashes),
                "raw_sha_authenticated_before_entry": True,
                "prior_receipt_decoded": False,
                "graph_decoded": False,
                "source_graph_imported": False,
                "source_graph_executed": False,
                "graph_traversed": False,
                "graph_edges_traversed": False,
                "actual_graph_evaluate_called": False,
                "actual_source_value_evaluated": False,
                "fixture_calls": 0,
                "integral_record_count": 0,
            },
        )
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_runner_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed_entry",
            "passed": False,
            "fixture": None,
            "integral_catalog": None,
            "traceback": traceback.format_exc(),
            "error": type(exc).__name__ + ": " + str(exc),
            "entry_written": False,
            "result_written": False,
        }
        print(json.dumps(output, ensure_ascii=False, sort_keys=True))
        return 1

    status = "returned"
    try:
        output = _consume(hashes, progress)
        output.update({"status": "returned", "passed": True, "traceback": None})
    except BaseException as exc:
        status = "failed"
        error_traceback = traceback.format_exc()
        progress["traceback"] = error_traceback
        output = {
            "status": status,
            "passed": False,
            "fixture": progress.get("fixture"),
            "fixture_raw": progress.get("fixture_raw"),
            "integral_catalog": progress.get("integral_catalog"),
            "integral_records_partial": progress.get("integral_records_partial", []),
            "traceback": error_traceback,
            "error": type(exc).__name__ + ": " + str(exc),
        }
    output.update(_base_result(source, hashes, source_sha, brief_sha, checkpoint_key, entry_sha, progress))
    try:
        result_sha = _write_exclusive(RESULT, output)
    except BaseException as exc:
        print(
            json.dumps(
                {
                    "component": COMPONENT,
                    "status": "failed_result_write",
                    "entry_sha256": entry_sha,
                    "error": type(exc).__name__ + ": " + str(exc),
                    "traceback": traceback.format_exc(),
                },
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 1

    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": status,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "source_checkpoint_sha256": source_sha,
                "pending_brief_sha256": brief_sha,
                "input_hashes": dict(hashes),
                "fixture_calls": int(progress.get("fixture_calls", 0)),
                "integral_record_count": len(progress.get("integral_records_partial", [])),
                "graph_decoded": bool(progress.get("graph_decoded")),
                "prior_receipt_decoded": bool(progress.get("prior_receipt_decoded")),
                "actual_graph_evaluate_called": False,
                "actual_source_value_evaluated": False,
                "graph_edges_traversed": False,
                "numeric_proof": False,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
