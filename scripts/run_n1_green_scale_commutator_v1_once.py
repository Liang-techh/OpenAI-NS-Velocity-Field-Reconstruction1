"""One-shot source-bound entry for the N1 green-scale commutator leaf.

The frozen manifest is authenticated by raw SHA before the exclusive entry is
created.  The parent cutoff report is deliberately not decoded during that
preflight.  After entry, the independent fixture is run once, then the new
leaf is called once with a mapping of verbatim source texts and the decoded
parent report.  None of the historical source producers are imported or
executed, and this runner does not claim a complete Kchi or modify a ledger.

The source manifest's input contract is intentionally explicit:

    fixture.run_fixture(leaf.audit_source_texts)
    leaf.audit_source_texts(
        source_texts["volterra"],
        source_texts["regular"],
        source_texts["first_raw"],
        parent=parent_report,
        source_labels=source_labels,
    )

``source_texts`` passed to the leaf contains the three actual Python source
texts under the role keys ``volterra``, ``regular``, and ``first_raw``;
the bridge document remains an authenticated provenance input.  The decoded
parent cutoff report is passed as ``parent_report``.
"""

from __future__ import annotations

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
COMPONENT = "n1_green_scale_commutator_v1"
PHASE = COMPONENT + "_pending"

BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_MODULE = "openai_ns_reconstruction.n1_green_scale_commutator_v1"
FIXTURE_MODULE = "openai_ns_reconstruction.n1_green_scale_commutator_fixture_v1"

LEAF_REL = "src/openai_ns_reconstruction/n1_green_scale_commutator_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_green_scale_commutator_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_green_scale_commutator_v1_once.py"
INNER_VOLTERRA_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
REGULAR_SERIES_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
)
ORIGINAL_RAW_REL = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"
PARENT_CUTOFF_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
)
BRIDGE_DOC_REL = "docs/CURRENT_N1_RAW_KCHI_RESIDUAL_BRIDGE_2026_11_02.md"

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    INNER_VOLTERRA_REL,
    REGULAR_SERIES_REL,
    ORIGINAL_RAW_REL,
    PARENT_CUTOFF_REL,
    BRIDGE_DOC_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
SOURCE_TEXT_INPUTS = (
    INNER_VOLTERRA_REL,
    REGULAR_SERIES_REL,
    ORIGINAL_RAW_REL,
)
PROVENANCE_INPUTS = (BRIDGE_DOC_REL,)
SOURCE_TEXT_ROLES = ("volterra", "regular", "first_raw")
LEAF_API = (
    "audit_source_texts(volterra_text, regular_text, first_raw_text, "
    "*, parent, source_labels)"
)
FIXTURE_API = "run_fixture(auditor)"
SHA256_HEX = frozenset("0123456789abcdefABCDEF")

# These three source pins were supplied by the parent handoff.  The other
# frozen inputs are pinned by the pending source manifest and reauthenticated
# against it before entry.
PINNED_INPUTS = {
    INNER_VOLTERRA_REL: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    REGULAR_SERIES_REL: "b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555",
    PARENT_CUTOFF_REL: "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856",
}

# Historical modules are forbidden even when a new leaf tries to load one
# under an alias.  The path-name check below complements the module-name
# check, while the empty namespace package prevents repository __init__.py
# execution.
FORBIDDEN_MODULE_BASES = frozenset(
    {
        "original_n1_first_raw_function",
        "n1_kchi_first_correction_v1",
        "n1_raw_kchi_residual_bridge_v1",
        "n1_raw_kchi_residual_bridge_fixture_v1",
        "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions",
        "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function",
    }
)
FORBIDDEN_FILE_NAMES = frozenset(name + ".py" for name in FORBIDDEN_MODULE_BASES)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


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


def _forbidden_module_names() -> list[str]:
    found: list[str] = []
    for name, module in sys.modules.items():
        base = name.rsplit(".", 1)[-1]
        module_file = getattr(module, "__file__", None)
        file_name = Path(module_file).name if isinstance(module_file, str) else ""
        if base in FORBIDDEN_MODULE_BASES or file_name in FORBIDDEN_FILE_NAMES:
            found.append(name)
    return sorted(set(found))


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    """Encode arbitrary leaf/fixture values without non-JSON objects."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "bytes", "bytes": len(value), "sha256": _sha_bytes(value)}

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
            items = [_json_safe(item, active) for item in value]
            return sorted(items, key=repr)
        return repr(value)
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
        ).encode("utf-8")
        + b"\n"
    )


def _write_exclusive(path: Path, value: Any) -> str:
    """Persist one JSON receipt with exclusive creation and fsync."""

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


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_bytes())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " path must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " path is missing or escapes repository: " + normalized,
    )
    return key, path


def _sha_text(value: Any, label: str) -> str:
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in SHA256_HEX for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate exactly the eight frozen files without decoding them."""

    entries: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            entries.append((item.get("path"), item.get("sha256"), f"source.files[{index}]"))
    else:
        for field in ("raw_file_sha256", "input_sha256", "input_hashes", "sha256"):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend(
                    (raw_path, expected, f"source.{field}")
                    for raw_path, expected in table.items()
                )
                break
    _require(entries, "source manifest contains no raw-file SHA entries")

    hashes: dict[str, str] = {}
    for raw_path, expected, label in entries:
        key, path = _repo_file(raw_path, label)
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        digest = _sha_text(expected, key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        actual = _sha_file(path)
        _require(actual.lower() == digest, "raw SHA mismatch: " + key)
        if key in PINNED_INPUTS:
            _require(actual.lower() == PINNED_INPUTS[key], "parent handoff pin changed: " + key)
        hashes[key] = actual

    _require(frozenset(hashes) == FROZEN_SET, "manifest must freeze exactly the eight required files")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint: Any = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        nested = brief.get(COMPONENT)
        if isinstance(nested, Mapping):
            checkpoint = nested.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        continuation = brief.get("continuation")
        if isinstance(continuation, Mapping):
            checkpoint = continuation.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")

    raw_path = checkpoint.get("path")
    _require(isinstance(raw_path, str), "pending brief source checkpoint path missing")
    path = Path(raw_path)
    if not path.is_absolute():
        path = ROOT / path
    _require(path.resolve() == SOURCE.resolve(), "pending brief points to another source")

    expected = checkpoint.get("sha256", checkpoint.get("manifest_sha256"))
    actual = _sha_file(SOURCE)
    _require(_sha_text(expected, "pending brief source checkpoint") == actual,
             "pending brief source checkpoint SHA mismatch")
    return actual


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


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(BRIEF.is_file(), "resume brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")

    # Only the brief and source manifest are decoded before entry.  Every
    # manifest-pinned content file, including the parent report, is hashed as
    # raw bytes here and is decoded only in _consume after entry.
    brief_raw = BRIEF.read_bytes()
    brief = _load_object(BRIEF, "resume brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")

    source = _load_object(SOURCE, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "source manifest phase mismatch")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")

    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "source manifest input contract missing")
    _require(
        contract.get("leaf_api") == LEAF_API,
        "leaf API contract mismatch",
    )
    _require(contract.get("fixture_api") == FIXTURE_API, "fixture API contract mismatch")
    contract_source_texts = tuple(contract.get("source_texts", ()))
    _require(
        contract_source_texts in (SOURCE_TEXT_ROLES, SOURCE_TEXT_INPUTS),
        "source-text input contract mismatch",
    )
    _require(contract.get("parent_report") == PARENT_CUTOFF_REL, "parent report contract mismatch")

    source_sha = _checkpoint_sha(brief)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha_bytes(brief_raw)


def _load_new_modules() -> tuple[Any, Any, list[str]]:
    """Load only the new modules under an empty namespace package."""

    preexisting_forbidden = _forbidden_module_names()
    _require(not preexisting_forbidden, "historical producer already loaded")

    source_root = str(SRC)
    package_name = "openai_ns_reconstruction"
    package_dir = SRC / package_name
    previous_package = sys.modules.get(package_name)
    previous_leaf = sys.modules.get(LEAF_MODULE)
    previous_fixture = sys.modules.get(FIXTURE_MODULE)

    namespace = ModuleType(package_name)
    namespace.__path__ = [str(package_dir)]
    namespace.__package__ = package_name
    namespace.__spec__ = importlib.machinery.ModuleSpec(
        package_name, loader=None, is_package=True
    )
    sys.modules[package_name] = namespace
    old_path = list(sys.path)
    if source_root not in sys.path:
        sys.path.insert(0, source_root)
    try:
        leaf_spec = importlib.util.spec_from_file_location(LEAF_MODULE, ROOT / LEAF_REL)
        fixture_spec = importlib.util.spec_from_file_location(FIXTURE_MODULE, ROOT / FIXTURE_REL)
        _require(leaf_spec is not None and leaf_spec.loader is not None,
                 "new leaf module spec is unavailable")
        _require(fixture_spec is not None and fixture_spec.loader is not None,
                 "new fixture module spec is unavailable")

        leaf = importlib.util.module_from_spec(leaf_spec)
        sys.modules[LEAF_MODULE] = leaf
        leaf_spec.loader.exec_module(leaf)

        fixture = importlib.util.module_from_spec(fixture_spec)
        sys.modules[FIXTURE_MODULE] = fixture
        fixture_spec.loader.exec_module(fixture)
    finally:
        sys.path[:] = old_path
        if previous_package is None:
            sys.modules.pop(package_name, None)
        else:
            sys.modules[package_name] = previous_package
        if previous_leaf is None:
            sys.modules.pop(LEAF_MODULE, None)
        else:
            sys.modules[LEAF_MODULE] = previous_leaf
        if previous_fixture is None:
            sys.modules.pop(FIXTURE_MODULE, None)
        else:
            sys.modules[FIXTURE_MODULE] = previous_fixture

    forbidden = _forbidden_module_names()
    _require(not forbidden, "new leaf/fixture imported a historical producer")
    return leaf, fixture, forbidden


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not a JSON object")
    result = dict(value)
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks), "fixture checks are missing")
    _require(all(item is True for item in checks.values()), "fixture has a failed check")
    for flag in (
        "historical_modules_called",
        "historical_producer_imported",
        "old_modules_called",
        "implementation_imported",
        "provider_executed",
        "numerical_residual_evaluated",
        "mpmath_used",
    ):
        if flag in result:
            _require(result[flag] is False, "fixture violated historical-source guard: " + flag)
    return _json_safe(result)


def _actual_scope(hashes: Mapping[str, str]) -> dict[str, Any]:
    return {
        "source_only": True,
        "source_text_inputs": [
            {"path": path, "sha256": hashes[path], "passed_verbatim": True}
            for path in SOURCE_TEXT_INPUTS
        ],
        "authenticated_provenance_inputs": [
            {"path": path, "sha256": hashes[path], "decoded": False}
            for path in PROVENANCE_INPUTS
        ],
        "parent_report_input": {
            "path": PARENT_CUTOFF_REL,
            "sha256": hashes[PARENT_CUTOFF_REL],
            "raw_sha_authenticated_before_entry": True,
            "decoded_after_entry_only": True,
        },
        "runner_rewrites_source_formula": False,
        "bare_G_replaced_by_parenthetical_G": False,
        "real_source_values_computed": False,
        "complete_kchi_proved": False,
        "ledger_modified": False,
        "full_goal_complete": False,
        "ancestor_math_calls": 0,
        "historical_producer_imported": bool(_forbidden_module_names()),
    }


def _base_result(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    entry_sha: str,
    calls: list[str],
    verification_invocations: int,
    actual_scope: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": f"{COMPONENT}_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "source_checkpoint_path": "tmp/" + COMPONENT + "/source.json",
        "source_checkpoint_sha256": source_sha,
        "source_manifest_sha256": source_sha,
        "brief_path": str(BRIEF),
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "raw_file_sha256": dict(hashes),
        "source_manifest": dict(source),
        "consumer_calls": list(calls),
        "verification_invocations": verification_invocations,
        "ancestor_math_calls": 0,
        "historical_producer_imported": bool(_forbidden_module_names()),
        "source_only": True,
        "actual_source_metadata_only": False,
        "complete_kchi_proved": False,
        "ledger_modified": False,
        "full_goal_complete": False,
        "actual_source_scope": dict(actual_scope),
        "runtime_model_metadata": None,
    }


def _consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    entry_sha: str,
    progress: dict[str, Any],
) -> dict[str, Any]:
    # Reauthenticate every frozen content byte after the exclusive entry.
    raw: dict[str, bytes] = {}
    for path in FROZEN_INPUTS:
        data = (ROOT / path).read_bytes()
        _require(_sha_bytes(data) == hashes[path], "input changed after entry: " + path)
        raw[path] = data

    leaf, fixture, forbidden = _load_new_modules()
    progress["historical_modules_imported"] = list(forbidden)

    # These are the first decodes of the actual inputs.  The runner passes the
    # source bytes verbatim as text and does not import either producer file.
    source_texts_by_path: dict[str, str] = {}
    for path in SOURCE_TEXT_INPUTS:
        try:
            source_texts_by_path[path] = raw[path].decode("utf-8")
        except UnicodeDecodeError as exc:
            raise RunnerInputError("source input is not UTF-8: " + path) from exc

    source_texts = {
        "volterra": source_texts_by_path[INNER_VOLTERRA_REL],
        "regular": source_texts_by_path[REGULAR_SERIES_REL],
        "first_raw": source_texts_by_path[ORIGINAL_RAW_REL],
    }

    try:
        parent_report_value = json.loads(raw[PARENT_CUTOFF_REL].decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError("parent cutoff report is not valid UTF-8 JSON") from exc
    _require(isinstance(parent_report_value, Mapping), "parent cutoff report is not a JSON object")
    parent_report = dict(parent_report_value)

    auditor = getattr(leaf, "audit_source_texts", None)
    _require(callable(auditor), "new leaf must export audit_source_texts(...)")
    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "fixture has no run_fixture(auditor, ...)")
    progress["verification_invocations"] = 1
    progress["consumer_calls"] = ["fixture.run_fixture(audit_source_texts)"]
    fixture_output = run_fixture(auditor)
    progress["fixture_raw"] = _json_safe(fixture_output)
    fixture_result = _validate_fixture(fixture_output)
    progress["fixture"] = fixture_result
    _require(
        not _forbidden_module_names(),
        "fixture imported a historical producer",
    )
    progress["fixture_inputs"] = {
        "fixture_owned_synthetic_texts": True,
        "fixture_source_texts_supplied_by_runner": False,
        "fixture_file_path": FIXTURE_REL,
        "fixture_file_sha256": hashes[FIXTURE_REL],
        "actual_source_texts_consumed": False,
    }

    progress["source_texts"] = {
        "roles": list(SOURCE_TEXT_ROLES),
        "paths": list(SOURCE_TEXT_INPUTS),
        "sha256": {path: hashes[path] for path in SOURCE_TEXT_INPUTS},
        "bytes": {path: len(raw[path]) for path in SOURCE_TEXT_INPUTS},
        "passed_verbatim": True,
        "decoded_after_entry": True,
    }
    progress["parent_report"] = {
        "path": PARENT_CUTOFF_REL,
        "sha256": hashes[PARENT_CUTOFF_REL],
        "bytes": len(raw[PARENT_CUTOFF_REL]),
        "decoded_after_entry": True,
    }

    source_labels = {
        "volterra": INNER_VOLTERRA_REL,
        "regular": REGULAR_SERIES_REL,
        "first_raw": ORIGINAL_RAW_REL,
    }
    progress["actual_leaf_calls"] = 1
    progress["consumer_calls"].append(
        "leaf.audit_source_texts(volterra_text, regular_text, first_raw_text, parent=..., source_labels=...)"
    )
    result = auditor(
        source_texts["volterra"],
        source_texts["regular"],
        source_texts["first_raw"],
        parent=parent_report,
        source_labels=source_labels,
    )
    progress["verification_invocations"] = 2
    encoded_result = _json_safe(result)
    # Preserve the complete actual leaf packet before validating its proof
    # flags, so a failed source proof is still written to the exclusive result.
    progress["actual_result"] = encoded_result
    _require(isinstance(result, Mapping), "leaf.audit_source_texts result must be a JSON object")
    _require(isinstance(encoded_result, Mapping), "encoded leaf result is not an object")
    _require(
        not _forbidden_module_names(),
        "new leaf imported a historical producer during audit",
    )
    _require(encoded_result.get("passed") is True, "actual source proof did not pass")
    checks = encoded_result.get("checks")
    _require(
        isinstance(checks, Mapping) and bool(checks),
        "actual source proof checks are missing",
    )
    _require(
        all(item is True for item in checks.values()),
        "actual source proof has a failed check",
    )

    # A source-bound bridge may report conditions and open premises, but this
    # runner refuses an accidental scope escalation in the leaf's own flags.
    for flag in (
        "complete_kchi_proved",
        "full_kchi_proved",
        "complete_Kchi_resolvent_provided",
        "full_Kchi_proved",
        "ledger_modified",
        "full_goal_complete",
    ):
        _require(encoded_result.get(flag) is not True,
                 "leaf attempted an out-of-scope claim: " + flag)

    actual_scope = _actual_scope(hashes)
    return {
        "status": "returned",
        "passed": True,
        "fixture": fixture_result,
        "leaf_result": encoded_result,
        "source_texts": progress["source_texts"],
        "fixture_inputs": progress["fixture_inputs"],
        "parent_report": progress["parent_report"],
        "actual_leaf_calls": 1,
        "ancestor_math_calls": 0,
        "historical_producer_imported": bool(_forbidden_module_names()),
        "source_only": True,
        "actual_source_metadata_only": False,
        "actual_source_scope": actual_scope,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "source_manifest_sha256": source_sha,
            "input_hashes": hashes,
            "raw_sha_authenticated_before_entry": True,
            "parent_report_decoded_before_entry": False,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "historical_producer_imported": False,
            "runtime_model_metadata": None,
        },
    )

    calls: list[str] = []
    verification_invocations = 0
    progress: dict[str, Any] = {
        "fixture": None,
        "actual_leaf_calls": 0,
        "consumer_calls": calls,
        "verification_invocations": verification_invocations,
    }
    actual_scope = _actual_scope(hashes)
    try:
        nested = _consume(source, hashes, source_sha, brief_sha, entry_sha, progress)
        verification_invocations = int(progress.get("verification_invocations", 0))
        calls = list(progress.get("consumer_calls", ()))
        output = _base_result(
            source, hashes, source_sha, brief_sha, entry_sha,
            calls, verification_invocations, actual_scope,
        )
        output.update(nested)
        output["consumer_calls"] = calls
        output["verification_invocations"] = verification_invocations
        output["error"] = None
        output["traceback"] = None
    except BaseException as exc:
        verification_invocations = int(progress.get("verification_invocations", 0))
        calls = list(progress.get("consumer_calls", ()))
        output = _base_result(
            source, hashes, source_sha, brief_sha, entry_sha,
            calls, verification_invocations, actual_scope,
        )
        output.update(
            {
                "status": "failed",
                "passed": False,
                "fixture": (
                    progress.get("fixture")
                    if progress.get("fixture") is not None
                    else progress.get("fixture_raw")
                ),
                "actual_result": progress.get("actual_result"),
                "progress": _json_safe(progress),
                "actual_leaf_calls": int(progress.get("actual_leaf_calls", 0)),
                "error": {"type": type(exc).__name__, "message": str(exc)},
                "traceback": traceback.format_exc(),
            }
        )

    result_sha = _write_exclusive(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": output["consumer_calls"],
                "verification_invocations": output["verification_invocations"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
