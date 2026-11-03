"""One-entry runner for the saved raw n=1 g5 C3 node records.

Raw SHA-256 authentication and the exclusive entry marker happen before any
new-module execution or receipt decoding.  Only the new leaf and fixture are
loaded; historical source modules are never imported.  The graph/check/binding
JSON is decoded only after the synthetic fixture passes, then the new leaf is
called exactly once.  The saved graph's C3 order is
``[value.node, Z.node, ZZ.node, ZZZ.node]``; the first two records are retained
as ordinary Z0/Z1 and the two charts remain separate.  No axis-owner identity
is claimed.
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
COMPONENT = "n1_saved_raw_g5_nodes_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_MODULE = PACKAGE + ".n1_saved_raw_g5_nodes_v1"
FIXTURE_MODULE = PACKAGE + ".n1_saved_raw_g5_nodes_fixture_v1"
LEAF_REL = "src/openai_ns_reconstruction/n1_saved_raw_g5_nodes_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_saved_raw_g5_nodes_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_saved_raw_g5_nodes_v1_once.py"
GRAPH_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
GRAPH_CHECK_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
GRAPH_SOURCE_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
CUTOFF_BINDING_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
C3_SOURCE_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions.py"
FROZEN_INPUTS = (LEAF_REL, FIXTURE_REL, RUNNER_REL, GRAPH_REL, GRAPH_CHECK_REL, GRAPH_SOURCE_REL, CUTOFF_BINDING_REL, C3_SOURCE_REL)
FROZEN_SET = frozenset(FROZEN_INPUTS)
GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
FIXTURE_SCHEMA = "n1_saved_raw_g5_nodes_fixture_v1"
LEAF_SCHEMA = "n1_saved_raw_g5_nodes_v1_result"
C3_ENCODING = ("value.node", "Z.node", "ZZ.node", "ZZZ.node")
CHARTS = ("core", "first")
SHA256_HEX = frozenset("0123456789abcdefABCDEF")
ALLOWED_MODULES = frozenset({PACKAGE, LEAF_MODULE, FIXTURE_MODULE})


class RunnerInputError(ValueError):
    pass


def _require(condition: bool, message: str) -> None:
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
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value), "sha256": _sha_bytes(value)}
    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        return repr(value)
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(_json_safe(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


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
    _require(isinstance(value, str) and len(value) == 64 and all(c in SHA256_HEX for c in value), "invalid SHA-256 for " + label)
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
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET, "manifest must contain exactly the eight frozen inputs")
    _require(hashes[GRAPH_REL] == GRAPH_SHA256, "saved raw graph SHA is not the accepted pin")
    return {key: hashes[key] for key in FROZEN_INPUTS}


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
    source = _decode_object(source_raw, "source manifest")
    _require(source.get("component") == COMPONENT and source.get("one_entry_only") is True, "source manifest guard failed")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "source manifest phase mismatch")
    return source, _manifest_hashes(source), source_sha, _sha_bytes(brief_raw)


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
    if "schema" in result:
        _require(result.get("schema") == FIXTURE_SCHEMA, "fixture schema mismatch")
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    checks = result.get("checks")
    if checks is not None:
        _require(isinstance(checks, Mapping) and bool(checks), "fixture checks must be a non-empty mapping")
        _require(all(item is True for item in checks.values()), "synthetic fixture contains a failed check")
    for key in ("historical_module_imported", "historical_modules_imported", "old_code_called", "old_receipt_decoded", "old_receipt_decompressed", "graph_decoded", "graph_decompressed", "leaf_called"):
        if key in result:
            _require(result[key] in (False, 0, None), "fixture reported forbidden activity: " + key)
    return _json_safe(result)


def _decode_graph(raw: bytes) -> dict[str, Any]:
    try:
        value = json.loads(gzip.decompress(raw))
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile, json.JSONDecodeError) as exc:
        raise RunnerInputError("authenticated raw graph is not valid gzip JSON") from exc
    _require(isinstance(value, Mapping), "saved raw graph is not a JSON object")
    return dict(value)


def _graph_g5_packets(graph: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    nodes = graph.get("exact_graph_nodes")
    bindings = graph.get("inner_Volterra_bindings")
    _require(isinstance(nodes, list) and nodes, "saved raw graph exact_graph_nodes missing")
    _require(isinstance(bindings, Mapping), "saved raw graph bindings missing")
    systems = bindings.get("systems")
    _require(isinstance(systems, Mapping) and set(systems) == set(CHARTS), "saved graph must expose core and first systems")
    packets: dict[str, dict[str, Any]] = {}
    for chart in CHARTS:
        system = systems[chart]
        _require(isinstance(system, Mapping), "saved chart system missing: " + chart)
        forcing = system.get("g")
        _require(isinstance(forcing, Mapping) and set(forcing) == {"4", "5", "6"}, "saved forcing rows changed: " + chart)
        encoding = forcing["5"]
        _require(isinstance(encoding, list) and len(encoding) == 4 and all(isinstance(item, int) and not isinstance(item, bool) for item in encoding), "g5 C3 encoding changed: " + chart)
        records: list[dict[str, Any]] = []
        for reference in encoding:
            _require(0 <= reference < len(nodes) and isinstance(nodes[reference], Mapping), "g5 node reference invalid: " + chart)
            records.append(dict(nodes[reference]))
        packets[chart] = {"encoding": list(encoding), "ordinary_Z0": encoding[0], "ordinary_Z1": encoding[1], "ordinary_Z_orders": [0, 1], "node_records": records, "C3_encoding": list(C3_ENCODING)}
    return packets


def _validate_graph_check(check: Mapping[str, Any]) -> None:
    _require(check.get("all_passed") is True, "saved graph check did not pass")
    _require(check.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is True, "saved graph check gate missing")


def _validate_binding(binding: Mapping[str, Any], graph: Mapping[str, Any]) -> None:
    graph_family = graph.get("source_family")
    binding_family = binding.get("source_family")
    if isinstance(graph_family, Mapping) and isinstance(binding_family, Mapping):
        _require(dict(binding_family) == dict(graph_family), "cutoff binding source family differs from graph")


def _validate_leaf_result(result: Any, expected: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    _require(isinstance(result, Mapping), "new leaf result is not a mapping")
    value = dict(result)
    _require(value.get("schema") == LEAF_SCHEMA, "new leaf result schema mismatch")
    _require(value.get("graph_hash") == GRAPH_SHA256, "new leaf graph hash changed")
    _require(value.get("file_path") == GRAPH_REL, "new leaf graph path changed")
    packets = value.get("g5")
    _require(isinstance(packets, Mapping) and set(packets) == set(CHARTS), "new leaf g5 packets must contain core and first")
    for chart in CHARTS:
        packet = packets[chart]
        _require(isinstance(packet, Mapping), "new leaf g5 packet missing: " + chart)
        _require(packet.get("saved_C3_node_ids") == expected[chart]["encoding"], "new leaf changed g5 C3 node IDs: " + chart)
        ordinary = packet.get("ordinary_Z")
        _require(isinstance(ordinary, Mapping) and set(ordinary) == {"Z0", "Z1"}, "new leaf lost ordinary g5 Z0/Z1: " + chart)
        for label, position in (("Z0", 0), ("Z1", 1)):
            item = ordinary[label]
            _require(isinstance(item, Mapping), "new leaf ordinary g5 packet is malformed: " + chart + "." + label)
            _require(item.get("ordinary_Z_label") == label and item.get("saved_C3_position") == position, "new leaf ordinary g5 position changed: " + chart + "." + label)
            _require(item.get("node_id") == expected[chart]["encoding"][position], "new leaf ordinary g5 node changed: " + chart + "." + label)
            _require(item.get("raw_record") == expected[chart]["node_records"][position], "new leaf did not preserve original g5 record: " + chart + "." + label)
    interpretation = value.get("interpretation")
    _require(isinstance(interpretation, Mapping), "new leaf interpretation is missing")
    _require(interpretation.get("ordinary_Z_labels") == ["Z0", "Z1"], "new leaf ordinary Z labels changed")
    _require(interpretation.get("references_are_symbolic_graph_nodes") is True, "new leaf node-reference marker missing")
    for key in ("point_values_computed", "live_owner_verified", "function_solution_handles_used", "axis_owner_claimed", "m15_claimed", "cancellation_claimed", "full_Kchi_claimed", "intervals_compared"):
        _require(interpretation.get(key) is False, "new leaf made an unsupported claim: " + key)
    return _json_safe(value)


def _consume(source: Mapping[str, Any], hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    after = _manifest_hashes(source)
    _require(dict(after) == dict(hashes), "manifest input changed after entry")
    raw = {path: (ROOT / path).read_bytes() for path in FROZEN_INPUTS}
    for path in FROZEN_INPUTS:
        _require(_sha_bytes(raw[path]) == hashes[path], "raw input changed after entry: " + path)
    progress["raw_sha_reauthenticated_after_entry"] = True
    try:
        leaf, fixture = _load_new_modules()
        progress["new_modules_imported"] = True
        run_fixture = getattr(fixture, "run_fixture", None)
        _require(callable(run_fixture), "fixture.run_fixture is missing or not callable")
        progress["fixture_calls"] = 1
        fixture_result = _validate_fixture(run_fixture(leaf))
        _require(not _forbidden_modules(), "historical module loaded during fixture")
        graph = _decode_graph(raw[GRAPH_REL])
        graph_check = _decode_object(raw[GRAPH_CHECK_REL], GRAPH_CHECK_REL)
        cutoff_binding = _decode_object(raw[CUTOFF_BINDING_REL], CUTOFF_BINDING_REL)
        _validate_graph_check(graph_check)
        _validate_binding(cutoff_binding, graph)
        expected = _graph_g5_packets(graph)
        graph_for_leaf = dict(graph)
        graph_for_leaf["graph_hash"] = hashes[GRAPH_REL]
        graph_for_leaf["file_path"] = GRAPH_REL
        consume = getattr(leaf, "consume_saved_raw_g5", None)
        _require(callable(consume), "new leaf consume_saved_raw_g5 API is missing or not callable")
        progress["leaf_calls"] = 1
        leaf_result = _validate_leaf_result(consume(graph_for_leaf), expected)
        return {"status": "returned", "passed": True, "fixture": fixture_result, "leaf_result": leaf_result, "graph_summary": {"source_family": _json_safe(graph.get("source_family")), "exact_graph_node_count": len(graph.get("exact_graph_nodes", [])), "charts": list(CHARTS), "C3_encoding": list(C3_ENCODING), "g5_original_nodes": _json_safe(expected)}, "graph_check_summary": {"all_passed": graph_check.get("all_passed"), "gate": graph_check.get("current_original_n1_true_inner_Volterra_function_adapter_installed")}, "cutoff_binding_summary": {"source_family": _json_safe(cutoff_binding.get("source_family")), "decoded_after_fixture": True}, "source_graph_imported": False, "source_graph_executed": False, "same_axis_owner_proved": False, "same_axis_owner_claimed": False, "same_source_owner_proved": False}
    finally:
        _cleanup_modules()


def _base_result(source: Mapping[str, Any], hashes: Mapping[str, str], source_sha: str, brief_sha: str, entry_sha: str, progress: Mapping[str, Any]) -> dict[str, Any]:
    return {"schema": COMPONENT + "_runner_result_v1", "component": COMPONENT, "phase": PHASE, "source_checkpoint_path": SOURCE.relative_to(ROOT).as_posix(), "source_checkpoint_sha256": source_sha, "source_manifest_sha256": source_sha, "brief_path": BRIEF.as_posix(), "brief_sha256": brief_sha, "entry_sha256": entry_sha, "input_hashes": dict(hashes), "raw_file_sha256": dict(hashes), "source_manifest": _json_safe(source), "one_entry_only": True, "do_not_reenter": True, "brief_modified": False, "fixture_calls": int(progress.get("fixture_calls", 0)), "leaf_calls": int(progress.get("leaf_calls", 0)), "progress": _json_safe(progress)}


def main() -> int:
    try:
        source, hashes, source_sha, brief_sha = _preflight()
    except BaseException as exc:
        print(json.dumps({"component": COMPONENT, "status": "preflight_failed", "error": str(exc), "traceback": traceback.format_exc()}, ensure_ascii=False, sort_keys=True))
        return 2
    try:
        entry_sha = _write_exclusive(ENTRY, {"schema": COMPONENT + "_entered_v1", "component": COMPONENT, "phase": PHASE, "entered": True, "one_entry_only": True, "do_not_reenter": True, "source_checkpoint_sha256": source_sha, "source_manifest_sha256": source_sha, "brief_sha256": brief_sha, "input_hashes": dict(hashes), "raw_sha_authenticated_before_entry": True, "source_graph_imported": False, "fixture_calls": 0, "leaf_calls": 0})
    except BaseException as exc:
        print(json.dumps({"component": COMPONENT, "status": "entry_failed", "error": str(exc), "traceback": traceback.format_exc()}, ensure_ascii=False, sort_keys=True))
        return 2
    progress: dict[str, Any] = {"fixture_calls": 0, "leaf_calls": 0}
    try:
        output = _base_result(source, hashes, source_sha, brief_sha, entry_sha, progress)
        output.update(_consume(source, hashes, progress))
        output.update({"error": None, "traceback": None})
        status = "returned"
    except BaseException as exc:
        output = _base_result(source, hashes, source_sha, brief_sha, entry_sha, progress)
        output.update({"status": "failed", "passed": False, "error": type(exc).__name__ + ": " + str(exc), "traceback": traceback.format_exc(), "source_graph_imported": False, "source_graph_executed": False, "same_axis_owner_proved": False, "same_axis_owner_claimed": False, "same_source_owner_proved": False})
        status = "failed"
    try:
        result_sha = _write_exclusive(RESULT, output)
    except BaseException as exc:
        print(json.dumps({"component": COMPONENT, "status": "result_write_failed", "entry_sha256": entry_sha, "error": str(exc), "traceback": traceback.format_exc()}, ensure_ascii=False, sort_keys=True))
        return 2
    print(json.dumps({"component": COMPONENT, "phase": PHASE, "status": status, "entry_sha256": entry_sha, "result_sha256": result_sha, "source_checkpoint_sha256": source_sha, "brief_sha256": brief_sha, "input_hashes": dict(hashes), "fixture_calls": int(progress.get("fixture_calls", 0)), "leaf_calls": int(progress.get("leaf_calls", 0))}, ensure_ascii=False, sort_keys=True))
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
