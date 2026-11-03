"""Run the v4 flat-sigma point/source-binding boundary exactly once.

This runner is a fresh v4 family.  It hashes every frozen input before the
exclusive entry, imports only the v4 leaf and fixture, runs the synthetic
fixture once, then authenticates and decodes only the prior g5 program receipt.
The v1, v2, and v3 failure archives are raw-only pins: they are hashed but never
read through gzip or JSON decoding.

The source AST evidence is reused from the saved v3 result; this runner does
not recompute the historical AST.  The saved Module-fragment provenance is
checked against pinned AST and raw source hashes; neither historical source
module is imported or executed.
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

COMPONENT = "n1_flat_sigma_point_kernel_v4"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_flat_sigma_point_kernel_v4.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_flat_sigma_point_kernel_fixture_v4.py"
RUNNER_REL = "scripts/run_n1_flat_sigma_point_kernel_v4_once.py"
PRIOR_RECEIPT_REL = "experiments/root_st073/n1_raw_g5_dag_program_v1_receipt.json.gz"
V1_FAILURE_REL = "experiments/root_st073/n1_flat_sigma_point_kernel_v1_failure.json.gz"
V2_FAILURE_REL = "experiments/root_st073/n1_flat_sigma_point_kernel_v2_failure.json.gz"
V3_FAILURE_REL = "experiments/root_st073/n1_flat_sigma_point_kernel_v3_failure.json.gz"
V3_RESULT_REL = "tmp/n1_flat_sigma_point_kernel_v3/result.json"
OUTER_SOURCE_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_outer_initial.py"
PULSE_SOURCE_REL = "experiments/root_st073/lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    PRIOR_RECEIPT_REL,
    V1_FAILURE_REL,
    V2_FAILURE_REL,
    V3_FAILURE_REL,
    V3_RESULT_REL,
    OUTER_SOURCE_REL,
    PULSE_SOURCE_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)

PRIOR_RECEIPT_SHA256 = "fa7108686ed6d6331bb8f09a12d86dba5ccd821d567aa45fbb1efc301f6f657f"
V1_FAILURE_SHA256 = "f746bdaa8430b264509cfac4ec9e2074efa19f0f9d1b6dd31479207396dd1a4f"
V2_FAILURE_SHA256 = "462d39067e91017f2d1b7c2616898be1c3fd785cbd41def61954c13169ec6874"
V3_FAILURE_SHA256 = "a56812e58b0c65fe979398706f3ef5948c5b84edb32539d2ac1ec9da913a14fa"
V3_RESULT_SHA256 = "b21f2de178bc17f74c4b68c041189bc61c0f5246eb6086294d88f8e30432bfea"
OUTER_SOURCE_SHA256 = "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657"
PULSE_SOURCE_SHA256 = "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7"
OUTER_AST_SHA256 = "4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c"
PULSE_AST_SHA256 = "a63c47d05049c86007b530f3095c515d724bef53092c5f7867d53e5bbaa06236"

PRIOR_RECEIPT_COMPONENT = "n1_raw_g5_dag_program_v1"
PRIOR_SCHEMA = "n1_raw_g5_dag_program_v1"
# This is a lineage label for the immediately preceding flat-sigma leaf.  It
# is not asserted from the archived g5 receipt's component field.
ARCHIVED_PREVIOUS_LABEL = "n1_flat_sigma_point_kernel_v3"
FIXTURE_SCHEMA = "n1_flat_sigma_point_kernel_fixture_v4"
LEAF_SCHEMA = "n1_flat_sigma_point_kernel_v4"
LEAF_MODULE = PACKAGE + ".n1_flat_sigma_point_kernel_v4"
FIXTURE_MODULE = PACKAGE + ".n1_flat_sigma_point_kernel_fixture_v4"
ALLOWED_MODULES = frozenset({PACKAGE, LEAF_MODULE, FIXTURE_MODULE})
SHA256_HEX = frozenset("0123456789abcdefABCDEF")

SIGMA_NODE_ID = 35949
SIGMA_PRIME_NODE_ID = 35953
ARGUMENT_NODE_ID = 35848


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
    """Record structures without calling arbitrary returned objects."""

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
        if key == PRIOR_RECEIPT_REL:
            _require(
                row.get("raw_sha256_before_decode") is True,
                "prior g5 receipt is not raw-pinned before decode",
            )
        if key in (V1_FAILURE_REL, V2_FAILURE_REL, V3_FAILURE_REL):
            _require(
                row.get("raw_sha256_before_decode") is True,
                "failure archive is not raw-pinned before decode: " + key,
            )
            _require(
                row.get("raw_sha256_only") is True,
                "failure archive is not marked raw-sha-only: " + key,
            )
        if key == V3_RESULT_REL:
            _require(
                row.get("raw_sha256_before_decode") is True,
                "v3 result is not raw-pinned before evidence read",
            )
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET, "manifest must contain exactly the ten frozen inputs")
    _require(hashes[PRIOR_RECEIPT_REL] == PRIOR_RECEIPT_SHA256, "prior receipt SHA is not the accepted pin")
    _require(hashes[V1_FAILURE_REL] == V1_FAILURE_SHA256, "v1 failure archive SHA is not the accepted pin")
    _require(hashes[V2_FAILURE_REL] == V2_FAILURE_SHA256, "v2 failure archive SHA is not the accepted pin")
    _require(hashes[V3_FAILURE_REL] == V3_FAILURE_SHA256, "v3 failure archive SHA is not the accepted pin")
    _require(hashes[V3_RESULT_REL] == V3_RESULT_SHA256, "v3 result SHA is not the accepted pin")
    _require(hashes[OUTER_SOURCE_REL] == OUTER_SOURCE_SHA256, "outer source SHA is not the accepted pin")
    _require(hashes[PULSE_SOURCE_REL] == PULSE_SOURCE_SHA256, "flat-pulse source SHA is not the accepted pin")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _validate_contract(source: Mapping[str, Any]) -> None:
    contract = source.get("input_contract")
    _require(contract is not None, "source input_contract is missing")
    _require(isinstance(contract, Mapping), "input_contract must be an object")
    _require(bool(contract), "source input_contract is empty")
    _require(
        contract.get("fixture_api") == "run_fixture(leaf)",
        "fixture API contract is missing or mismatched",
    )
    _require(
        contract.get("binding_api")
        == "bind_flat_sigma_point_kernel(prior_program, source_digests)",
        "binding API contract is missing or mismatched",
    )
    _require(
        contract.get("prior_receipt_path") == PRIOR_RECEIPT_REL,
        "prior receipt path is missing or mismatched",
    )
    _require(
        contract.get("prior_program_path") == ["result", "program"],
        "prior program path is missing or mismatched",
    )
    _require(contract.get("v1_failure_path") == V1_FAILURE_REL, "v1 failure path mismatch")
    _require(contract.get("v2_failure_path") == V2_FAILURE_REL, "v2 failure path mismatch")
    _require(contract.get("v3_failure_path") == V3_FAILURE_REL, "v3 failure path mismatch")
    _require(contract.get("v3_result_path") == V3_RESULT_REL, "v3 result path mismatch")
    _require(contract.get("v1_failure_raw_only") is True, "v1 failure raw-only contract is disabled")
    _require(contract.get("v2_failure_raw_only") is True, "v2 failure raw-only contract is disabled")
    _require(contract.get("v3_failure_raw_only") is True, "v3 failure raw-only contract is disabled")
    _require(contract.get("v3_result_raw_before_decode") is True, "v3 result raw-before-decode marker is disabled")
    _require(
        contract.get("source_ast_algorithm")
        == "source_fragment_ast_parse_module_dump",
        "source AST algorithm marker is missing or mismatched",
    )
    _require(
        contract.get("authenticate_raw_gzip_before_decode") is True,
        "raw gzip authentication marker is missing or disabled",
    )
    _require(
        contract.get("no_graph_import") is True,
        "no_graph_import marker is missing or disabled",
    )
    _require(
        contract.get("source_ast_evidence")
        == "reuse_v3_result_progress_source_fragment_provenance",
        "saved v3 AST evidence reuse marker is missing or mismatched",
    )


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


def _brief_checkpoint(brief: Mapping[str, Any], record: Mapping[str, Any]) -> Mapping[str, Any]:
    candidates: list[Any] = [
        record.get("source_checkpoint"),
        brief.get("source_checkpoint"),
        brief.get(COMPONENT + "_source_checkpoint"),
    ]
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        candidates.append(nested.get("source_checkpoint"))
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping):
        candidates.append(continuation.get("source_checkpoint"))
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            return candidate
    raise RunnerInputError("pending brief source_checkpoint is missing")


def _validate_pending_brief(source_sha: str) -> tuple[str, str]:
    _require(BRIEF.is_file(), "pending resume brief is missing")
    brief_raw = BRIEF.read_bytes()
    brief_sha = _sha_bytes(brief_raw)
    brief = _decode_object(brief_raw, "pending resume brief")
    record = _record_for(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending brief one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending brief re-entry guard missing")
    checkpoint = _brief_checkpoint(brief, record)
    raw_path = checkpoint.get("path")
    _require(isinstance(raw_path, str) and raw_path, "pending brief source checkpoint path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "pending brief source checkpoint path must be repository-relative",
    )
    checkpoint_key = relative.as_posix()
    expected_key = SOURCE.relative_to(ROOT).as_posix()
    _require(checkpoint_key == expected_key, "pending brief points to another source checkpoint")
    _require((ROOT / relative).resolve() == SOURCE.resolve(), "pending brief checkpoint resolves incorrectly")
    _require(
        _sha_text(checkpoint.get("sha256"), "pending source checkpoint") == source_sha,
        "pending brief source checkpoint SHA mismatch",
    )
    return brief_sha, checkpoint_key


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


def _revalidate_after_entry(
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    checkpoint_key: str,
    progress: dict[str, Any],
) -> None:
    """Re-authenticate every manifest row after entry and before imports."""

    source_raw = SOURCE.read_bytes()
    _require(_sha_bytes(source_raw) == source_sha, "source manifest changed after entry")
    current_source = _decode_object(source_raw, "source manifest after entry")
    _require(current_source.get("component") == COMPONENT, "source manifest component changed after entry")
    _require(current_source.get("phase") == PHASE, "source manifest phase changed after entry")
    _require(current_source.get("one_entry_only") is True, "source manifest one-entry guard changed after entry")
    _require(current_source.get("do_not_reenter") is True, "source manifest re-entry guard changed after entry")
    _validate_contract(current_source)
    current_hashes = _manifest_hashes(current_source)
    _require(dict(current_hashes) == dict(hashes), "manifest file hashes changed after entry")
    current_brief_sha, current_checkpoint_key = _validate_pending_brief(source_sha)
    _require(current_brief_sha == brief_sha, "pending brief changed after entry")
    _require(current_checkpoint_key == checkpoint_key, "pending brief checkpoint changed after entry")
    progress["post_entry_manifest_rehash_verified"] = True
    progress["post_entry_checkpoint_verified"] = True


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
    def checked(
        name: str,
        globals_: Any = None,
        locals_: Any = None,
        fromlist: tuple[str, ...] | None = None,
        level: int = 0,
    ):
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
    _require(receipt.get("component") == PRIOR_RECEIPT_COMPONENT, "prior receipt component mismatch")
    result = receipt.get("result")
    _require(isinstance(result, Mapping), "prior receipt result is missing")
    _require(result.get("status") == "returned", "prior receipt result is not returned")
    _require(result.get("passed") is True, "prior receipt result did not pass")
    program = result.get("program")
    _require(isinstance(program, Mapping), "prior receipt result.program is missing")
    _require(program.get("schema") == PRIOR_SCHEMA, "prior program schema mismatch")
    _require(
        program.get("status") == "compiled" or program.get("compile_status") == "compiled",
        "prior program is not explicitly compiled",
    )
    return program


def _saved_v3_source_evidence(
    saved_result: Mapping[str, Any],
    hashes: Mapping[str, str],
) -> tuple[dict[str, dict[str, str]], dict[str, Any]]:
    """Reuse v3's saved Module-fragment evidence without recomputing AST."""

    _require(saved_result.get("component") == "n1_flat_sigma_point_kernel_v3", "saved v3 result component mismatch")
    _require(saved_result.get("status") == "failed", "saved v3 result is not failed")
    _require(saved_result.get("passed") is False, "saved v3 result passed unexpectedly")
    _require(saved_result.get("source_fragment_ast_verified") is True, "saved v3 AST evidence is not verified")
    _require(saved_result.get("source_identity_verified") is False, "saved v3 result claimed source identity")
    _require(saved_result.get("numeric_proof") is False, "saved v3 result claimed numeric proof")
    _require(saved_result.get("actual_source_value_evaluated") is False, "saved v3 result evaluated source values")
    _require(saved_result.get("actual_point_values_computed") is False, "saved v3 result computed actual point values")

    progress = saved_result.get("progress")
    _require(isinstance(progress, Mapping), "saved v3 progress is missing")
    _require(progress.get("source_fragment_ast_verified") is True, "saved v3 progress AST evidence is not verified")
    _require(progress.get("source_fragment_ast_recomputed") in (False, None), "saved v3 claimed AST recomputation")
    provenance = progress.get("source_fragment_provenance")
    _require(isinstance(provenance, Mapping), "saved v3 source fragment provenance is missing")

    saved_hashes = saved_result.get("input_hashes")
    _require(isinstance(saved_hashes, Mapping), "saved v3 input hashes are missing")
    for path in (PRIOR_RECEIPT_REL, OUTER_SOURCE_REL, PULSE_SOURCE_REL):
        _require(saved_hashes.get(path) == hashes[path], "saved v3 input hash mismatch: " + path)

    saved_manifest = saved_result.get("source_manifest")
    _require(isinstance(saved_manifest, Mapping), "saved v3 source manifest is missing")
    _require(saved_manifest.get("component") == "n1_flat_sigma_point_kernel_v3", "saved v3 manifest component mismatch")
    saved_files = saved_manifest.get("files")
    _require(isinstance(saved_files, list), "saved v3 manifest files are missing")
    saved_file_hashes = {
        row.get("path"): row.get("sha256")
        for row in saved_files
        if isinstance(row, Mapping)
    }
    for path in (PRIOR_RECEIPT_REL, OUTER_SOURCE_REL, PULSE_SOURCE_REL):
        _require(saved_file_hashes.get(path) == hashes[path], "saved v3 manifest source hash mismatch: " + path)

    expected = {
        "sigma": ("stable_sigma", OUTER_AST_SHA256),
        "sigma_prime": ("sigma_jets", PULSE_AST_SHA256),
    }
    digests: dict[str, dict[str, str]] = {}
    for role, (function_name, ast_sha) in expected.items():
        evidence = provenance.get(role)
        _require(isinstance(evidence, Mapping), "saved v3 AST evidence is missing: " + role)
        _require(evidence.get("algorithm") == "source_fragment_ast_parse_module_dump", "saved v3 AST algorithm mismatch: " + role)
        _require(evidence.get("ast_container") == "Module", "saved v3 AST container mismatch: " + role)
        _require(evidence.get("ast_dump_defaults") is True, "saved v3 AST dump settings mismatch: " + role)
        _require(evidence.get("fragment_semantics_match") is True, "saved v3 fragment semantics mismatch: " + role)
        _require(evidence.get("function") == function_name, "saved v3 function mismatch: " + role)
        _require(evidence.get("AST_sha256") == ast_sha, "saved v3 AST SHA mismatch: " + role)
        if role == "sigma":
            module = Path(OUTER_SOURCE_REL).name
            module_sha = hashes[OUTER_SOURCE_REL]
        else:
            module = Path(PULSE_SOURCE_REL).name
            module_sha = hashes[PULSE_SOURCE_REL]
        digests[role] = {
            "module": module,
            "function": function_name,
            "module_sha256": module_sha,
            "AST_sha256": ast_sha,
        }

    return digests, {
        "sigma": dict(provenance["sigma"]),
        "sigma_prime": dict(provenance["sigma_prime"]),
        "source_fragment_ast_verified": True,
        "source_fragment_ast_recomputed": False,
        "source_ast_evidence_reused_from": V3_RESULT_REL,
        "saved_v3_result_status": saved_result.get("status"),
    }


def _validate_binding(value: Any, supplied_digests: Mapping[str, Any]) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "binding result is not a mapping")
    result = dict(value)
    _require(result.get("schema") == LEAF_SCHEMA, "binding schema mismatch")
    _require(result.get("status") == "bound", "flat-sigma v4 binding did not return bound")

    _require(result.get("source_identity_verified") is False, "binding claimed source identity verification")
    _require(
        result.get("caller_supplied_source_records_match") is True,
        "caller-supplied source records were not matched",
    )
    _require(result.get("prior_program_executable") is False, "prior program was marked executable")
    _require(result.get("prior_program_execution_supported") is False, "prior execution support was claimed")
    _require(result.get("program_executable") is False, "bound program was marked executable")
    _require(result.get("compiled_program_executable") is False, "compiled program execution was claimed")
    _require(result.get("point_evaluation_only") is True, "point-only scope was not retained")
    _require(result.get("interval") is False, "binding incorrectly claimed interval output")
    _require(result.get("certified") is False, "binding incorrectly claimed certification")
    _require(
        result.get("numeric_equivalence_certified") is False,
        "binding incorrectly claimed numeric equivalence",
    )
    _require(result.get("full_Navier_Stokes") is False, "binding incorrectly claimed full Navier-Stokes output")
    claims = result.get("claims")
    _require(isinstance(claims, Mapping), "binding claims record is missing")
    _require(claims.get("certified") is False, "binding claims certified output")
    _require(claims.get("interval") is False, "binding claims interval output")
    _require(claims.get("numeric_equivalence_certified") is False, "binding claims numeric equivalence")
    _require(claims.get("ast_recomputed") is False, "binding claims a fresh AST recomputation")
    _require(claims.get("source_identity_verified") is False, "binding claims source identity verification")
    _require(claims.get("prior_program_executable") is False, "binding claims prior execution")
    identities = result.get("source_identities")
    _require(isinstance(identities, Mapping), "binding source identity records are missing")
    for role, node_id in (("sigma", SIGMA_NODE_ID), ("sigma_prime", SIGMA_PRIME_NODE_ID)):
        supplied = supplied_digests.get(role)
        _require(isinstance(supplied, Mapping), f"caller {role} source digest is missing")
        identity = identities.get(role)
        _require(isinstance(identity, Mapping), f"binding {role} source identity is missing")
        for key in ("module", "module_sha256", "AST_sha256"):
            _require(identity.get(key) == supplied.get(key), f"binding {role} source {key} mismatch")
        if "function" in supplied:
            _require(identity.get("function") == supplied.get("function"), f"binding {role} source function mismatch")
    _require(result.get("bound_argument_node_id") == ARGUMENT_NODE_ID, "bound argument node changed")
    bound_nodes = result.get("bound_nodes")
    _require(isinstance(bound_nodes, Mapping), "bound node records are missing")
    for role, node_id in (("sigma", SIGMA_NODE_ID), ("sigma_prime", SIGMA_PRIME_NODE_ID)):
        row = bound_nodes.get(role)
        _require(isinstance(row, Mapping), f"bound {role} record is missing")
        _require(row.get("node_id") == node_id, f"bound {role} node id changed")
        _require(row.get("argument") == ARGUMENT_NODE_ID, f"bound {role} argument changed")
    return _json_safe(result)


def _scope_summary(progress: Mapping[str, Any]) -> dict[str, Any]:
    """Describe the narrow provenance scope without implying AST or graph proof."""

    binding_success = bool(progress.get("binding_success"))
    evidence_reused = bool(progress.get("source_ast_evidence_reused"))
    return {
        "target_source_node_ids": [SIGMA_NODE_ID, SIGMA_PRIME_NODE_ID],
        "target_source_node_count": 2,
        "shared_argument_node_required": ARGUMENT_NODE_ID,
        "only_target_source_nodes_located": True,
        "record_resolution": (
            "target wrappers are preferred; bare records are an allowed fallback; "
            "conflicting wrapper, record, and alias metadata is rejected"
        ),
        "saved_v3_source_ast_summary_reused": evidence_reused,
        "saved_ast_assertion_reused_and_target_records_match": (
            binding_success and evidence_reused
        ),
        "node_source_hash_ast_provenance_match": False,
        "current_ast_recomputed": False,
        "numeric_equivalence_certified": False,
        "full_graph_binding": False,
        "scope_description": (
            "Only the two target source nodes are located against the pinned prior "
            "records, with the shared argument node required. The AST fields are "
            "summaries reused from the saved v3 result; this is not a current AST "
            "recomputation, numeric-equivalence proof, or full-graph binding."
        ),
    }


def _flags(progress: Mapping[str, Any]) -> dict[str, Any]:
    binding_success = bool(progress.get("binding_success"))
    fragment_verified = bool(progress.get("source_fragment_ast_verified"))
    evidence_reused = bool(progress.get("source_ast_evidence_reused"))
    return {
        "raw_inputs_hashed": bool(progress.get("raw_inputs_hashed")),
        "entry_written": bool(progress.get("entry_written")),
        "post_entry_manifest_rehash_verified": bool(progress.get("post_entry_manifest_rehash_verified")),
        "post_entry_checkpoint_verified": bool(progress.get("post_entry_checkpoint_verified")),
        "new_modules_imported": bool(progress.get("new_modules_imported")),
        "fixture_calls": int(progress.get("fixture_calls", 0)),
        "fixture_passed": bool(progress.get("fixture_passed")),
        "prior_receipt_decoded": bool(progress.get("prior_receipt_decoded")),
        "prior_receipt_raw_hash_rechecked": bool(progress.get("prior_receipt_raw_hash_rechecked")),
        "v3_result_raw_hash_rechecked": bool(progress.get("v3_result_raw_hash_rechecked")),
        "v3_result_loaded": bool(progress.get("v3_result_loaded")),
        "source_ast_evidence_reused": bool(progress.get("source_ast_evidence_reused")),
        "source_fragment_ast_recomputed": False,
        "binding_calls": int(progress.get("binding_calls", 0)),
        "binding_success": binding_success,
        "source_fragment_ast_verified": fragment_verified,
        "caller_supplied_source_records_match": bool(progress.get("caller_supplied_source_records_match")),
        "source_identity_verified": bool(progress.get("source_identity_verified")),
        "saved_ast_assertion_reused_and_target_records_match": bool(
            binding_success and evidence_reused
        ),
        "node_source_hash_ast_provenance_match": False,
        "numeric_proof": False,
        "actual_source_value_evaluated": False,
        "actual_point_values_computed": False,
        "specified_Z_evaluated": False,
        "v1_failure_raw_sha256_only": bool(progress.get("v1_failure_raw_sha256_only")),
        "v2_failure_raw_sha256_only": bool(progress.get("v2_failure_raw_sha256_only")),
        "v3_failure_raw_sha256_only": bool(progress.get("v3_failure_raw_sha256_only")),
        "v1_failure_decompressed": False,
        "v2_failure_decompressed": False,
        "v3_failure_decompressed": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "graph_traversed": False,
        "old_producer_imported": False,
        "graph_import_used": False,
    }


def _consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    leaf, fixture = _load_new_modules()
    progress["new_modules_imported"] = True
    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "fixture.run_fixture is missing or not callable")
    progress["fixture_calls"] = 1
    raw_fixture = run_fixture(leaf)
    progress["fixture_raw"] = _json_safe(raw_fixture)
    fixture_result = _validate_fixture(raw_fixture)
    progress["fixture"] = fixture_result
    progress["fixture_passed"] = True
    _require(not _forbidden_modules(), "historical module loaded during fixture")

    # Reuse the completed v3 AST evidence after the fresh fixture, without
    # parsing source text or recomputing either AST digest.
    v3_result_path = ROOT / V3_RESULT_REL
    raw_v3_result = v3_result_path.read_bytes()
    _require(_sha_bytes(raw_v3_result) == hashes[V3_RESULT_REL], "v3 result changed after preflight")
    progress["v3_result_raw_hash_rechecked"] = True
    saved_v3_result = _decode_object(raw_v3_result, V3_RESULT_REL)
    progress["v3_result_loaded"] = True
    digests, provenance = _saved_v3_source_evidence(saved_v3_result, hashes)
    progress["source_fragment_provenance"] = provenance
    progress["source_fragment_ast_verified"] = True
    progress["source_ast_evidence_reused"] = True

    # Decode the pinned prior receipt only after the fresh fixture and saved
    # v3 evidence have passed their checks.
    receipt_path = ROOT / PRIOR_RECEIPT_REL
    raw_receipt = receipt_path.read_bytes()
    _require(_sha_bytes(raw_receipt) == hashes[PRIOR_RECEIPT_REL], "prior receipt changed after preflight")
    progress["prior_receipt_raw_hash_rechecked"] = True
    receipt = _decode_gzip_json(raw_receipt, PRIOR_RECEIPT_REL)
    progress["prior_receipt_decoded"] = True
    program = _prior_program(receipt)

    bind = getattr(leaf, "bind_flat_sigma_point_kernel", None)
    _require(callable(bind), "v4 leaf binding API is missing or not callable")
    progress["binding_calls"] = 1
    raw_binding = bind(prior_program=program, source_digests=digests)
    progress["binding_raw"] = _json_safe(raw_binding)
    binding = _validate_binding(raw_binding, digests)
    progress["binding"] = binding
    progress["binding_success"] = True
    progress["caller_supplied_source_records_match"] = True
    _require(not _forbidden_modules(), "historical module loaded by new files")

    flags = _flags(progress)
    return {
        "status": "returned",
        "passed": True,
        "fixture": fixture_result,
        "binding": binding,
        "prior_receipt_component": receipt.get("component"),
        "archived_previous_label": ARCHIVED_PREVIOUS_LABEL,
        "prior_program_schema": program.get("schema"),
        "prior_program_path": ["result", "program"],
        "caller_supplied_source_digests": digests,
        "source_fragment_provenance": provenance,
        "source_binding_scope": "pinned prior node records plus source module hashes and reused AST evidence only",
        "source_binding_is_metadata_only": True,
        "source_ast_evidence_reused": True,
        "source_fragment_ast_recomputed": False,
        "numeric_authentication": False,
        "actual_source_value_evaluated": False,
        "actual_point_values_computed": False,
        "v1_failure_archive": {
            "path": V1_FAILURE_REL,
            "sha256": hashes[V1_FAILURE_REL],
            "raw_sha256_only": True,
            "decompressed": False,
        },
        "v2_failure_archive": {
            "path": V2_FAILURE_REL,
            "sha256": hashes[V2_FAILURE_REL],
            "raw_sha256_only": True,
            "decompressed": False,
        },
        "v3_failure_archive": {
            "path": V3_FAILURE_REL,
            "sha256": hashes[V3_FAILURE_REL],
            "raw_sha256_only": True,
            "decompressed": False,
        },
        "v3_result_evidence_archive": {
            "path": V3_RESULT_REL,
            "sha256": hashes[V3_RESULT_REL],
            "raw_sha256_before_decode": True,
            "ast_evidence_reused": True,
            "ast_recomputed": False,
        },
        **flags,
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
        "archived_previous_label": ARCHIVED_PREVIOUS_LABEL,
        "source_checkpoint_path": SOURCE.relative_to(ROOT).as_posix(),
        "pending_brief_path": BRIEF.relative_to(ROOT.parent.parent).as_posix(),
        "pending_brief_sha256": brief_sha,
        "pending_brief_checkpoint_path": checkpoint_key,
        "source_checkpoint_sha256": source_sha,
        "input_hashes": dict(hashes),
        "one_entry_only": True,
        "do_not_reenter": True,
        "brief_modified": False,
        "raw_gzip_authenticated_before_decode": bool(progress.get("prior_receipt_raw_hash_rechecked")),
        "failure_archives_raw_sha256_only": True,
        "failure_archives_decompressed": False,
        "source_ast_evidence_reused": bool(progress.get("source_ast_evidence_reused")),
        "source_fragment_ast_recomputed": False,
        "v3_result_evidence_loaded": bool(progress.get("v3_result_loaded")),
        "summary": _scope_summary(progress),
        "source_manifest": _json_safe(source),
        "entry_sha256": entry_sha,
        **_flags(progress),
        "progress": _json_safe(progress),
    }


def main() -> int:
    try:
        source, hashes, source_sha, brief_sha, checkpoint_key = _preflight()
    except BaseException as exc:
        preflight_progress = {
            "raw_inputs_hashed": False,
            "entry_written": False,
            "post_entry_manifest_rehash_verified": False,
            "post_entry_checkpoint_verified": False,
            "new_modules_imported": False,
            "fixture_calls": 0,
            "fixture_passed": False,
            "prior_receipt_decoded": False,
            "prior_receipt_raw_hash_rechecked": False,
            "v3_result_raw_hash_rechecked": False,
            "v3_result_loaded": False,
            "source_ast_evidence_reused": False,
            "source_fragment_ast_recomputed": False,
            "binding_calls": 0,
            "binding_success": False,
            "source_fragment_ast_verified": False,

            "caller_supplied_source_records_match": False,
            "source_identity_verified": False,
            "v1_failure_raw_sha256_only": False,
            "v2_failure_raw_sha256_only": False,
            "v3_failure_raw_sha256_only": False,
            "v1_failure_decompressed": False,
            "v2_failure_decompressed": False,
            "v3_failure_decompressed": False,
        }
        print(
            json.dumps(
                {
                    "schema": COMPONENT + "_runner_result_v1",
                    "component": COMPONENT,
                    "phase": PHASE,
                    "archived_previous_label": ARCHIVED_PREVIOUS_LABEL,
                    "status": "failed_preflight",
                    "passed": False,
                    "error": type(exc).__name__ + ": " + str(exc),
                    "traceback": traceback.format_exc(),
                    **_flags(preflight_progress),
                },
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 1

    progress: dict[str, Any] = {
        "raw_inputs_hashed": True,
        "entry_written": False,
        "post_entry_manifest_rehash_verified": False,
        "post_entry_checkpoint_verified": False,
        "new_modules_imported": False,
        "fixture_calls": 0,
        "fixture_passed": False,
        "prior_receipt_decoded": False,
        "prior_receipt_raw_hash_rechecked": False,
        "v3_result_raw_hash_rechecked": False,
        "v3_result_loaded": False,
        "source_ast_evidence_reused": False,
        "source_fragment_ast_recomputed": False,
        "binding_calls": 0,
        "binding_success": False,
        "source_fragment_ast_verified": False,
        "caller_supplied_source_records_match": False,
        "source_identity_verified": False,
        "v1_failure_raw_sha256_only": True,
        "v2_failure_raw_sha256_only": True,
        "v3_failure_raw_sha256_only": True,
        "v1_failure_decompressed": False,
        "v2_failure_decompressed": False,
        "v3_failure_decompressed": False,
    }
    try:
        entry_sha = _write_exclusive(
            ENTRY,
            {
                "schema": COMPONENT + "_entered_v1",
                "component": COMPONENT,
                "phase": PHASE,
                "archived_previous_label": ARCHIVED_PREVIOUS_LABEL,
                "entered": True,
                "one_entry_only": True,
                "do_not_reenter": True,
                "source_checkpoint_sha256": source_sha,
                "pending_brief_sha256": brief_sha,
                "pending_brief_checkpoint_path": checkpoint_key,
                "input_hashes": dict(hashes),
                "raw_sha_authenticated_before_entry": True,
                **_flags(progress),
            },
        )
        progress["entry_written"] = True
    except BaseException as exc:
        print(
            json.dumps(
                {
                    "schema": COMPONENT + "_runner_result_v1",
                    "component": COMPONENT,
                    "phase": PHASE,
                    "archived_previous_label": ARCHIVED_PREVIOUS_LABEL,
                    "status": "failed_entry",
                    "passed": False,
                    "error": type(exc).__name__ + ": " + str(exc),
                    "traceback": traceback.format_exc(),
                    **_flags(progress),
                },
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 1

    status = "returned"
    try:
        _revalidate_after_entry(hashes, source_sha, brief_sha, checkpoint_key, progress)
        output = _consume(hashes, progress)
        output.update({"status": "returned", "passed": True})
    except BaseException as exc:
        status = "failed"
        output = {
            "status": status,
            "passed": False,
            "error": type(exc).__name__ + ": " + str(exc),
            "traceback": traceback.format_exc(),
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
                **_flags(progress),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())




