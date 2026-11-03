"""Run the fresh raw-g5 DAG compiler exactly once.

This runner is intentionally a provenance and symbolic-program boundary.  It
authenticates every manifest input before importing the new leaf or fixture,
writes the exclusive entry record, runs ``fixture.run_fixture(leaf)`` once,
and only then decodes the authenticated saved graph and g5 extraction receipt.
The real graph is passed to ``leaf.compile_program`` once.  The runner never
evaluates a graph node, materializes a point value, or claims source identity,
physical m15 meaning, or a complete Kchi resolvent.

The two prior gzip files are raw-byte pins.  Their contents are opaque until
the preflight and exclusive-entry gates have passed.  Only the new leaf and
fixture are loaded in an empty ``openai_ns_reconstruction`` namespace; old
package modules and historical source modules are rejected.
"""

from __future__ import annotations

import builtins
import gzip
import hashlib
import importlib.machinery
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
SRC = ROOT / "src"
PACKAGE = "openai_ns_reconstruction"
PACKAGE_DIR = SRC / PACKAGE

COMPONENT = "n1_raw_g5_dag_program_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_raw_g5_dag_program_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_raw_g5_dag_program_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_raw_g5_dag_program_v1_once.py"
GRAPH_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
G5_RECEIPT_REL = "experiments/root_st073/n1_saved_raw_g5_nodes_v1_receipt.json.gz"
FROZEN_INPUTS = (LEAF_REL, FIXTURE_REL, RUNNER_REL, GRAPH_REL, G5_RECEIPT_REL)
FROZEN_SET = frozenset(FROZEN_INPUTS)
RECEIPTS = frozenset({GRAPH_REL, G5_RECEIPT_REL})

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
G5_RECEIPT_SHA256 = "3eb58012ba9057403991bbb035d08e1bed2224c38840cb38ebe31ea76d9e4f69"
PIN = {GRAPH_REL: GRAPH_SHA256, G5_RECEIPT_REL: G5_RECEIPT_SHA256}

GRAPH_NODE_COUNT = 49024
CHARTS = ("core", "first")
C3_ENCODING = ("value.node", "Z.node", "ZZ.node", "ZZZ.node")
FIXTURE_SCHEMA = "n1_raw_g5_dag_program_fixture_v1"
EXPECTED_G5 = {
    "core": (37863, 38521, 38687, 39097),
    "first": (40991, 42163, 42379, 42787),
}

LEAF_MODULE = PACKAGE + ".n1_raw_g5_dag_program_v1"
FIXTURE_MODULE = PACKAGE + ".n1_raw_g5_dag_program_fixture_v1"
ALLOWED_MODULES = frozenset({PACKAGE, LEAF_MODULE, FIXTURE_MODULE})
SHA256_HEX = frozenset("0123456789abcdefABCDEF")


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
    """Preserve returned structures without evaluating source data."""

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
    _require(not relative.is_absolute() and ".." not in relative.parts, label + " path must be repository-relative")
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(), label + " path is missing or escapes repository: " + key)
    return key, path


def _checkpoint_path(value: Any) -> Path:
    _require(isinstance(value, str) and value, "source checkpoint path missing")
    path = Path(value)
    if path.is_absolute():
        path = path.resolve()
        _require(path == SOURCE.resolve(), "pending brief points to another source checkpoint")
        return path
    return _repo_file(value, "source checkpoint")[1]


def _record_for(brief: Mapping[str, Any]) -> Mapping[str, Any]:
    if brief.get("component") == COMPONENT:
        return brief
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        return nested
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping) and continuation.get("component") == COMPONENT:
        return continuation
    return brief


def _brief_checkpoint(brief: Mapping[str, Any]) -> Mapping[str, Any]:
    candidates: list[Any] = [brief.get("source_checkpoint"), brief.get(COMPONENT + "_source_checkpoint")]
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        candidates.append(nested.get("source_checkpoint"))
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping):
        candidates.append(continuation.get("source_checkpoint"))
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            return candidate
    raise RunnerInputError("pending source checkpoint is missing")


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
        if key in PIN:
            _require(actual == PIN[key], "pinned input changed: " + key)
            _require(row.get("raw_sha256_before_decode") is True, "gzip input is not marked raw-pinned before decode: " + key)
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET, "manifest must contain exactly the five frozen inputs")
    _require(hashes[GRAPH_REL] == GRAPH_SHA256, "saved graph SHA is not the accepted pin")
    _require(hashes[G5_RECEIPT_REL] == G5_RECEIPT_SHA256, "g5 receipt SHA is not the accepted pin")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _validate_contract(source: Mapping[str, Any]) -> None:
    """Validate optional manifest API declarations without inventing them."""

    contract = source.get("input_contract")
    if contract is None:
        return
    _require(isinstance(contract, Mapping), "input_contract must be an object")
    if "leaf_api" in contract:
        _require(contract.get("leaf_api") == "compile_program(graph, roots, graph_hash)", "leaf API contract mismatch")
    if "fixture_api" in contract:
        _require(contract.get("fixture_api") == "run_fixture(leaf)", "fixture API contract mismatch")
    if "graph_path" in contract:
        _require(contract.get("graph_path") == GRAPH_REL, "graph contract path mismatch")
    if "g5_receipt_path" in contract:
        _require(contract.get("g5_receipt_path") == G5_RECEIPT_REL, "g5 receipt contract path mismatch")
    if "graph_sha256" in contract:
        _require(contract.get("graph_sha256") == GRAPH_SHA256, "graph contract hash mismatch")
    if "g5_receipt_sha256" in contract:
        _require(contract.get("g5_receipt_sha256") == G5_RECEIPT_SHA256, "g5 receipt contract hash mismatch")
    if "authenticate_raw_gzip_before_decode" in contract:
        _require(contract.get("authenticate_raw_gzip_before_decode") is True, "raw gzip authentication contract is disabled")


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered or has a result")
    _require(BRIEF.is_file() and SOURCE.is_file(), "pending brief or source manifest is missing")

    brief_raw = BRIEF.read_bytes()
    brief = _decode_object(brief_raw, "resume brief")
    record = _record_for(brief)
    _require(record.get("component") == COMPONENT and record.get("phase") == PHASE, "pending brief component/phase mismatch")
    _require(record.get("one_entry_only") is True and record.get("do_not_reenter") is True, "pending one-entry guard missing")
    _require(record.get("goal_active") is True and record.get("full_goal_complete") is False, "pending goal guard failed")

    checkpoint = _brief_checkpoint(brief)
    source_raw = SOURCE.read_bytes()
    source_sha = _sha_bytes(source_raw)
    _require(_checkpoint_path(checkpoint.get("path")) == SOURCE.resolve(), "source checkpoint path mismatch")
    _require(_sha_text(checkpoint.get("sha256"), "source checkpoint") == source_sha, "source checkpoint SHA mismatch")
    preimage_path = checkpoint.get("preimage_path")
    preimage_sha = checkpoint.get("preimage_sha256")
    if preimage_path is not None or preimage_sha is not None:
        _require(isinstance(preimage_path, str) and isinstance(preimage_sha, str), "incomplete source preimage checkpoint")
        _require(_sha_text(preimage_sha, "source preimage") == _sha_file(Path(preimage_path)), "source preimage SHA mismatch")

    source = _decode_object(source_raw, "source manifest")
    _require(source.get("component") == COMPONENT and source.get("one_entry_only") is True, "source manifest guard failed")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "source manifest phase mismatch")
    _validate_contract(source)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha_bytes(brief_raw)


def _forbidden_modules() -> list[str]:
    found: list[str] = []
    package_root = PACKAGE_DIR.resolve()
    experiments_root = (ROOT / "experiments").resolve()
    for name, module in list(sys.modules.items()):
        if name in ALLOWED_MODULES:
            continue
        if name == PACKAGE or name.startswith(PACKAGE + "."):
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


def _guarded_import(original: Any):
    def checked(name: str, globals_: Any = None, locals_: Any = None, fromlist: tuple[str, ...] | None = None, level: int = 0):
        resolved = name
        if level:
            package = globals_.get("__package__", "") if isinstance(globals_, Mapping) else ""
            try:
                resolved = importlib.util.resolve_name(name, package)
            except (ImportError, AttributeError):
                resolved = name
        if resolved == PACKAGE:
            for item in fromlist or ():
                if item != "*" and PACKAGE + "." + item not in ALLOWED_MODULES:
                    raise RunnerInputError("historical package import blocked: " + item)
        elif resolved.startswith(PACKAGE + ".") and resolved not in ALLOWED_MODULES:
            raise RunnerInputError("historical package import blocked: " + resolved)
        return original(name, globals_, locals_, fromlist, level)

    return checked


def _cleanup_modules() -> None:
    for name in list(sys.modules):
        if name == PACKAGE or name.startswith(PACKAGE + "."):
            sys.modules.pop(name, None)


def _load_new_modules() -> tuple[ModuleType, ModuleType]:
    _require(not _forbidden_modules(), "historical package or source module already loaded")
    namespace = ModuleType(PACKAGE)
    namespace.__path__ = [str(PACKAGE_DIR)]
    namespace.__package__ = PACKAGE
    namespace.__spec__ = importlib.machinery.ModuleSpec(PACKAGE, loader=None, is_package=True)
    sys.modules[PACKAGE] = namespace

    old_path = list(sys.path)
    original_import = builtins.__import__
    builtins.__import__ = _guarded_import(original_import)
    try:
        if str(SRC) not in sys.path:
            sys.path.insert(0, str(SRC))
        loaded: list[ModuleType] = []
        for module_name, relative in ((LEAF_MODULE, LEAF_REL), (FIXTURE_MODULE, FIXTURE_REL)):
            path = ROOT / relative
            spec = importlib.util.spec_from_file_location(module_name, path)
            _require(spec is not None and spec.loader is not None, "new module spec is unavailable: " + relative)
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
            loaded.append(module)
        _require(not _forbidden_modules(), "historical module loaded by new files")
        return loaded[0], loaded[1]
    except BaseException:
        _cleanup_modules()
        raise
    finally:
        builtins.__import__ = original_import
        sys.path[:] = old_path


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not a mapping")
    result = dict(value)
    _require(result.get("schema") == FIXTURE_SCHEMA, "synthetic fixture schema mismatch")
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    checks = result.get("checks")
    if checks is not None:
        _require(isinstance(checks, Mapping) and bool(checks), "fixture checks must be a non-empty mapping")
        _require(all(item is True for item in checks.values()), "synthetic fixture contains a failed check")
    for key in (
        "historical_module_imported",
        "historical_modules_imported",
        "old_code_called",
        "old_receipt_decoded",
        "old_receipt_decompressed",
        "graph_decoded",
        "graph_decompressed",
        "leaf_called",
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


def _saved_g5_result(receipt: Mapping[str, Any]) -> dict[str, Any]:
    """Read the prior runner's direct ``receipt['result']`` payload.

    The accepted archive stores ``leaf_result`` directly under ``result``;
    there is no intermediate ``result['output']`` field.
    """

    result = receipt.get("result")
    _require(isinstance(result, Mapping), "g5 extraction receipt result is missing")
    _require(result.get("status") == "returned", "g5 extraction receipt is not returned")
    return dict(result)


def _validate_saved_g5(receipt: Mapping[str, Any]) -> dict[str, Any]:
    result = _saved_g5_result(receipt)
    _require(result.get("component") == "n1_saved_raw_g5_nodes_v1", "unexpected g5 extraction component")
    leaf_result = result.get("leaf_result")
    _require(isinstance(leaf_result, Mapping), "g5 extraction leaf result is missing")
    _require(leaf_result.get("schema") == "n1_saved_raw_g5_nodes_v1_result", "unexpected g5 extraction schema")
    _require(leaf_result.get("graph_hash") == GRAPH_SHA256, "g5 extraction graph hash mismatch")
    _require(leaf_result.get("graph_node_count") == GRAPH_NODE_COUNT, "g5 extraction graph node count mismatch")
    interpretation = leaf_result.get("interpretation")
    _require(isinstance(interpretation, Mapping), "g5 extraction interpretation is missing")
    for key in (
        "point_values_computed",
        "live_owner_verified",
        "function_solution_handles_used",
        "axis_owner_claimed",
        "m15_claimed",
        "cancellation_claimed",
        "full_Kchi_claimed",
        "intervals_compared",
    ):
        _require(interpretation.get(key) is False, "g5 extraction made unsupported claim: " + key)
    graph_summary = result.get("graph_summary")
    _require(isinstance(graph_summary, Mapping), "g5 extraction graph summary is missing")
    _require(graph_summary.get("C3_encoding") == list(C3_ENCODING), "g5 extraction C3 encoding labels changed")
    g5 = leaf_result.get("g5")
    _require(isinstance(g5, Mapping) and set(g5) == set(CHARTS), "g5 extraction charts are incomplete")
    return {"result": result, "leaf_result": dict(leaf_result)}


def _decode_graph(raw: bytes) -> dict[str, Any]:
    graph = _decode_gzip_json(raw, GRAPH_REL)
    nodes = graph.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) == GRAPH_NODE_COUNT, "saved graph node count changed")
    if "graph_hash" in graph:
        _require(graph.get("graph_hash") == GRAPH_SHA256, "saved graph embedded hash mismatch")
    return graph


def _extract_roots(saved: Mapping[str, Any]) -> dict[str, Any]:
    """Extract only the four authenticated ordinary-Z root handles.

    ``handle`` and ``node_id`` are the same saved graph reference.  The raw
    node record is copied for provenance but is never evaluated; no child
    node is followed here.
    """

    leaf_result = saved["leaf_result"]
    saved_g5 = leaf_result["g5"]
    roots: dict[str, Any] = {}
    for chart in CHARTS:
        expected = EXPECTED_G5[chart]
        packet = saved_g5[chart]
        _require(isinstance(packet, Mapping), "g5 extraction packet is malformed: " + chart)
        _require(tuple(packet.get("saved_C3_node_ids") or ()) == expected, "g5 C3 IDs changed: " + chart)
        ordinary = packet.get("ordinary_Z")
        _require(isinstance(ordinary, Mapping) and set(ordinary) == {"Z0", "Z1"}, "ordinary roots missing: " + chart)

        chart_roots: dict[str, Any] = {}
        for label, position in (("Z0", 0), ("Z1", 1)):
            item = ordinary[label]
            _require(isinstance(item, Mapping), "ordinary root is not a mapping: " + chart + "." + label)
            node_id = item.get("node_id")
            _require(type(node_id) is int and node_id == expected[position], "ordinary root node ID changed: " + chart + "." + label)
            node_path = item.get("node_path")
            _require(node_path == ["exact_graph_nodes", node_id], "ordinary root path changed: " + chart + "." + label)
            source_path = item.get("source_path")
            expected_source_path = ["inner_Volterra_bindings", "systems", chart, "g", "5"]
            _require(source_path == expected_source_path, "ordinary root source path changed: " + chart + "." + label)
            _require(item.get("ordinary_Z_label") == label, "ordinary root label changed: " + chart + "." + label)
            _require(item.get("saved_C3_position") == position, "ordinary root C3 position changed: " + chart + "." + label)
            _require("raw_record" in item, "ordinary root raw record missing: " + chart + "." + label)
            raw_record = item["raw_record"]
            chart_roots[label] = {
                "chart": chart,
                "ordinary_Z": label,
                "ordinary_Z_order": position,
                "handle": node_id,
                "node_id": node_id,
                "node_path": list(node_path),
                "raw_record": _json_safe(raw_record),
                "source_path": list(source_path),
                "graph_hash": GRAPH_SHA256,
            }
        roots[chart] = chart_roots
    return roots


def _validate_program_scope(value: Any) -> None:
    """Reject explicit positive claims outside this symbolic compiler scope."""

    # The new leaf is expected to return its mapping program, but an opaque
    # caller-owned program object must remain recordable rather than being
    # rejected merely for being opaque.  Scope checks below apply when the
    # returned value exposes named claim fields.
    if not isinstance(value, Mapping):
        return
    forbidden_true = (
        "eval_executed",
        "graph_executed",
        "numeric_evaluation",
        "numeric_values_computed",
        "point_values_computed",
        "source_values_computed",
        "same_source_owner_proved",
        "same_axis_owner_proved",
        "m15_claimed",
        "full_Kchi_claimed",
        "authenticated",
    )
    for key in forbidden_true:
        if key in value:
            _require(value.get(key) is not True, "compile_program made unsupported claim: " + key)


def _consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    progress: dict[str, Any],
) -> dict[str, Any]:
    source_raw_after = SOURCE.read_bytes()
    _require(_sha_bytes(source_raw_after) == source_sha, "source manifest changed after entry")
    source_after = _decode_object(source_raw_after, "source manifest after entry")
    _validate_contract(source_after)
    _require(_manifest_hashes(source_after) == dict(hashes), "manifest input changed after entry")

    raw: dict[str, bytes] = {}
    for path in FROZEN_INPUTS:
        data = (ROOT / path).read_bytes()
        _require(_sha_bytes(data) == hashes[path], "raw input changed after entry: " + path)
        raw[path] = data
    progress["raw_sha_reauthenticated_after_entry"] = True

    try:
        leaf, fixture = _load_new_modules()
        progress["new_modules_imported"] = True
        run_fixture = getattr(fixture, "run_fixture", None)
        _require(callable(run_fixture), "fixture.run_fixture is missing or not callable")
        progress["fixture_calls"] = 1
        fixture_result = _validate_fixture(run_fixture(leaf))
        progress["fixture"] = fixture_result
        _require(not _forbidden_modules(), "historical module loaded during fixture")

        graph = _decode_graph(raw[GRAPH_REL])
        saved_receipt = _decode_gzip_json(raw[G5_RECEIPT_REL], G5_RECEIPT_REL)
        saved = _validate_saved_g5(saved_receipt)
        roots = _extract_roots(saved)
        root_order = (("core", "Z0"), ("core", "Z1"), ("first", "Z0"), ("first", "Z1"))
        root_ids = [roots[chart][label]["node_id"] for chart, label in root_order]
        progress["decoded_after_fixture"] = True
        progress["graph_summary"] = {
            "graph_hash": GRAPH_SHA256,
            "graph_node_count": len(graph["exact_graph_nodes"]),
            "charts": list(CHARTS),
            "C3_encoding": list(C3_ENCODING),
            "g5_C3_node_ids": {chart: list(EXPECTED_G5[chart]) for chart in CHARTS},
        }
        receipt_result = saved_receipt.get("result")
        progress["g5_receipt"] = {
            "path": G5_RECEIPT_REL,
            "sha256": hashes[G5_RECEIPT_REL],
            "component": saved_receipt.get("component"),
            "status": receipt_result.get("status") if isinstance(receipt_result, Mapping) else None,
        }
        progress["root_handles"] = roots
        progress["root_ids"] = root_ids
        progress["root_order"] = [chart + "." + label for chart, label in root_order]

        compile_program = getattr(leaf, "compile_program", None)
        _require(callable(compile_program), "new leaf.compile_program API is missing or not callable")
        progress["compile_program_calls"] = 1
        program = compile_program(graph, root_ids, GRAPH_SHA256)
        # Save the complete return before scope validation so an unsupported
        # claim cannot erase the returned symbolic program boundary.
        progress["program_return"] = _json_safe(program)
        _validate_program_scope(program)
        if isinstance(program, Mapping):
            if "roots" in program:
                _require(program.get("roots") == root_ids, "compile_program changed root ID order")
            if "graph_hash" in program:
                _require(program.get("graph_hash") == GRAPH_SHA256, "compile_program changed graph hash")
        return {
            "status": "returned",
            "passed": True,
            "fixture": fixture_result,
            "graph_summary": progress["graph_summary"],
            "g5_receipt": progress["g5_receipt"],
            "root_handles": roots,
            "root_ids": root_ids,
            "root_order": [chart + "." + label for chart, label in root_order],
            "program": _json_safe(program),
            "compile_program_calls": 1,
            "fixture_calls": 1,
            "source_graph_evaluated": False,
            "numeric_values_computed": False,
            "same_source_claimed": False,
            "m15_claimed": False,
            "full_Kchi_claimed": False,
            "unknown_and_opaque_boundaries_retained": True,
        }
    finally:
        _cleanup_modules()


def _base_result(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    entry_sha: str,
    progress: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": COMPONENT + "_runner_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "source_checkpoint_path": SOURCE.relative_to(ROOT).as_posix(),
        "source_checkpoint_sha256": source_sha,
        "source_manifest_sha256": source_sha,
        "brief_path": BRIEF.as_posix(),
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "raw_file_sha256": dict(hashes),
        "source_manifest": _json_safe(source),
        "one_entry_only": True,
        "do_not_reenter": True,
        "brief_modified": False,
        "raw_gzip_authenticated_before_decode": True,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "ancestor_called": False,
        "fixture_calls": int(progress.get("fixture_calls", 0)),
        "compile_program_calls": int(progress.get("compile_program_calls", 0)),
        "progress": _json_safe(progress),
    }


def main() -> int:
    try:
        source, hashes, source_sha, brief_sha = _preflight()
    except BaseException as exc:
        print(
            json.dumps(
                {"component": COMPONENT, "status": "preflight_failed", "error": str(exc), "traceback": traceback.format_exc()},
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 2

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
                "source_manifest_sha256": source_sha,
                "brief_sha256": brief_sha,
                "input_hashes": dict(hashes),
                "raw_gzip_authenticated_before_decode": True,
                "source_graph_imported": False,
                "fixture_calls": 0,
                "compile_program_calls": 0,
                "ancestor_called": False,
            },
        )
    except BaseException as exc:
        print(
            json.dumps(
                {"component": COMPONENT, "status": "entry_failed", "error": str(exc), "traceback": traceback.format_exc()},
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 2

    progress: dict[str, Any] = {"fixture_calls": 0, "compile_program_calls": 0}
    try:
        output = _base_result(source, hashes, source_sha, brief_sha, entry_sha, progress)
        output.update(_consume(source, hashes, source_sha, progress))
        output.update({"error": None, "traceback": None})
        status = "returned"
    except BaseException as exc:
        output = _base_result(source, hashes, source_sha, brief_sha, entry_sha, progress)
        output.update(
            {
                "status": "failed",
                "passed": False,
                "error": type(exc).__name__ + ": " + str(exc),
                "traceback": traceback.format_exc(),
                "source_graph_imported": False,
                "source_graph_executed": False,
                "ancestor_called": False,
                "same_source_claimed": False,
                "m15_claimed": False,
                "full_Kchi_claimed": False,
            }
        )
        status = "failed"

    try:
        result_sha = _write_exclusive(RESULT, output)
    except BaseException as exc:
        print(
            json.dumps(
                {"component": COMPONENT, "status": "result_write_failed", "entry_sha256": entry_sha, "error": str(exc), "traceback": traceback.format_exc()},
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 2

    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": status,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "source_checkpoint_sha256": source_sha,
                "brief_sha256": brief_sha,
                "input_hashes": dict(hashes),
                "fixture_calls": int(progress.get("fixture_calls", 0)),
                "compile_program_calls": int(progress.get("compile_program_calls", 0)),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
