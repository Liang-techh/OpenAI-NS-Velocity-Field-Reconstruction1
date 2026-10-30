"""Run the bounded F-gap index-3 source mapping capture exactly once.

This runner is the persistence boundary for
``f_gap_z3_row_mapping_capture_v1``.  It authenticates exactly five raw
repository inputs before writing the exclusive entry receipt:

* the new capture leaf;
* this runner;
* the previously published frame main/gap proof result;
* the authenticated Rp gzip archive; and
* the authenticated Rp source-check receipt.

The runner does not parse those three prior source/result inputs.  The leaf
authenticates them again before parsing and decodes the gzip only after its raw
SHA-256 matches.  There is no synthetic fixture call.  The only post-entry
leaf call is the real bounded index-3 source consumer.  Missing or protected
mapping evidence remains unresolved; no zero or other value is synthesized.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_z3_row_mapping_capture_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER_REL = "scripts/run_f_gap_z3_row_mapping_capture_once.py"
LEAF_REL = "src/openai_ns_reconstruction/f_gap_z3_row_mapping_capture.py"
PRIOR_REL = "tmp/frame_main_gap_source_proof_v1/result.json"
FRAME_ARCHIVE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
RECEIPT_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source_check.json"
)
SOURCE_INPUT_KEY = (
    "lei_ren_part1_paper_compliant_current_original_"
    "n1_pulse_pressure_C3_functions.json.gz"
)

REQUIRED_INPUTS = frozenset(
    {LEAF_REL, RUNNER_REL, PRIOR_REL, FRAME_ARCHIVE_REL, RECEIPT_REL}
)
FIXED_SHA256 = {
    PRIOR_REL: "af33465e04cd8901060227f4ab75f072efdd5c4d2fd59893a7ff9add27404144",
    FRAME_ARCHIVE_REL: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    RECEIPT_REL: "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "source_values_selected",
        "source_value_selected",
        "source_values_computed",
        "numeric_values_computed",
        "numeric_difference_evaluated",
        "numeric_error_computed",
        "numeric_evaluation",
        "admissible_for_numeric_use",
        "graph_accessed",
        "graph_traversed",
        "graph_edges_followed",
        "exact_graph_nodes_accessed",
        "children_traversed",
        "child_rows_indexed",
        "source_constructors_called",
        "source_checker_called",
        "old_helper_imported",
        "old_constructor_called",
        "ancestor_replayed",
        "ordinary_Z3_semantics_proved",
        "m3_ordinary_Z3_semantics_proved",
        "ordinary_Z_transport_proved",
        "F_source_identity_proved",
        "physical_radius_proved",
        "global_closure_proved",
        "full_goal_complete",
        "error_removed",
        "old_error_removed",
        "old_errors_removed",
        "ledger_212_213_183_changed",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return value.as_posix()
    return repr(value)


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(_json_safe(value))
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_bytes())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _repo_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    candidate = (ROOT / relative).resolve()
    _require(
        candidate.is_relative_to(ROOT.resolve()) and candidate.is_file(),
        label + " is missing or escapes repository: " + normalized,
    )
    return key, candidate


def _manifest_entries(source: Mapping[str, Any]) -> list[tuple[Any, Any, str]]:
    entries: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(
                isinstance(item, Mapping),
                f"source.files[{index}] is not an object",
            )
            entries.append(
                (item.get("path"), item.get("sha256"), f"source.files[{index}]")
            )
    else:
        for field in (
            "sha256",
            "input_hashes",
            "input_sha256",
            "raw_input_sha256",
            "code_sha256",
        ):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend(
                    (raw_path, expected, f"source.{field}")
                    for raw_path, expected in table.items()
                )
    _require(entries, "source manifest contains no direct raw SHA entries")
    return entries


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    entries = _manifest_entries(source)
    _require(
        len(entries) == len(REQUIRED_INPUTS),
        "source manifest must contain exactly five raw inputs",
    )
    hashes: dict[str, str] = {}
    for raw_path, expected, label in entries:
        key, path = _repo_path(raw_path, label)
        _require(key in REQUIRED_INPUTS, "manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        digest = expected.lower()
        actual = _sha256_file(path)
        _require(actual.lower() == digest, "raw SHA mismatch: " + key)
        fixed = FIXED_SHA256.get(key)
        if fixed is not None:
            _require(digest == fixed, "fixed SHA pin changed: " + key)
        hashes[key] = actual
    _require(
        frozenset(hashes) == REQUIRED_INPUTS,
        "manifest must freeze exactly the five required raw inputs",
    )
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    expected = checkpoint.get("sha256")
    _require(isinstance(raw_path, str), "source checkpoint path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(), "pending brief points elsewhere")
    actual = _sha256_file(SOURCE)
    _require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and actual.lower() == expected.lower(),
        "pending brief source checkpoint SHA-256 mismatch",
    )
    return actual


def _guard_false_scope(value: Any, label: str) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                raise RunnerInputError(
                    f"{label} made forbidden positive claim at {key_text}"
                )
            _guard_false_scope(item, label)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _guard_false_scope(item, label)


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(
        not ENTRY.exists() and not RESULT.exists(),
        "one-entry component already entered",
    )
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = _load_object(BRIEF, "pending brief")

    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief active/no-reentry guards changed",
    )
    one_shot = brief.get("one_entry_only")
    if one_shot is None:
        one_shot = brief.get("one_shot")
    _require(one_shot is True, "pending brief one-shot guard missing")

    # Authenticate source.json bytes before parsing the manifest or writing the
    # exclusive entered receipt.
    source_sha = _checkpoint_sha(brief)
    source = _load_object(SOURCE, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    if "schema" in source:
        _require(
            source.get("schema") == COMPONENT + "_checkpoint_v1",
            "source manifest schema mismatch",
        )
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")
    _require(
        source.get("synthetic_only") is False
        or source.get("non_synthetic") is True
        or source.get("real_source_only") is True,
        "source manifest must declare a non-synthetic source consumer",
    )
    for key in (
        "no_ancestor_replay",
        "no_old_math",
        "no_old_producer",
        "no_graph_traversal",
        "no_source_constructor",
        "no_source_checker",
    ):
        if key in source:
            _require(source[key] is True, "source manifest guard changed: " + key)
    _guard_false_scope(source, "source manifest")
    hashes = _manifest_hashes(source)
    brief_sha = hashlib.sha256(brief_raw).hexdigest()
    return source, hashes, source_sha, brief_sha


def _import_leaf() -> Any:
    leaf = (ROOT / LEAF_REL).resolve()
    _require(leaf.is_file(), "new z3-row mapping leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_z3_row_mapping_capture_v1_new_leaf", leaf
    )
    _require(spec is not None and spec.loader is not None, "new leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _leaf_callable(module: Any) -> tuple[str, Callable[..., Any]]:
    for name in ("build_f_gap_z3_row_mapping_capture", "recover"):
        function = getattr(module, name, None)
        if callable(function):
            return name, function
    raise RunnerInputError(
        "new leaf must expose build_f_gap_z3_row_mapping_capture or recover"
    )


def _validate_result(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "new leaf did not return an object")
    result = dict(value)
    _guard_false_scope(result, "new leaf result")
    _require(result.get("schema") == COMPONENT, "new leaf schema mismatch")
    _require(result.get("index") == 3, "new leaf did not retain index 3")
    _require(
        result.get("admissible_for_numeric_use") is False,
        "new leaf admitted numeric use",
    )
    _require(
        result.get("source_values_selected") is False,
        "new leaf selected source values",
    )
    flags = result.get("flags")
    if isinstance(flags, Mapping):
        for key in (
            "m3_ordinary_Z3_semantics_proved",
            "ordinary_Z3_semantics_proved",
            "ordinary_Z_transport_proved",
            "F_source_identity_proved",
            "physical_radius_proved",
            "numeric_values_computed",
            "numeric_difference_evaluated",
            "numeric_error_computed",
            "global_closure_proved",
            "full_goal_complete",
        ):
            if key in flags:
                _require(flags[key] is False, "new leaf flag must remain false: " + key)
    return result


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "manifest_paths_exact": sorted(hashes),
            "new_leaf_imported": False,
            "leaf_api_invocations": 0,
            "synthetic_fixture_invocations": 0,
            "source_consumer_invocations": 0,
            "gzip_decoded_by_runner": False,
            "graph_accessed": False,
            "source_values_selected": False,
            "ordinary_Z3_semantics_proved": False,
            "physical_radius_proved": False,
            "numeric_error_computed": False,
            "global_closure_proved": False,
            "runtime_model_metadata": None,
        },
    )

    leaf_imported = False
    call_name: str | None = None
    source_result: Any = None
    try:
        module = _import_leaf()
        leaf_imported = True
        call_name, consumer = _leaf_callable(module)
        source_result = consumer(
            ROOT,
            snapshot_path=PRIOR_REL,
            snapshot_sha256=FIXED_SHA256[PRIOR_REL],
            frame_archive_path=FRAME_ARCHIVE_REL,
            frame_archive_sha256=FIXED_SHA256[FRAME_ARCHIVE_REL],
            receipt_path=RECEIPT_REL,
            receipt_sha256=FIXED_SHA256[RECEIPT_REL],
            source_input_key=SOURCE_INPUT_KEY,
        )
        source_result = _validate_result(source_result)
        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "manifest_paths_exact": sorted(hashes),
            "leaf_imported": leaf_imported,
            "leaf_callable": call_name,
            "leaf_api_invocations": 1,
            "synthetic_fixture_invocations": 0,
            "source_consumer_invocations": 1,
            "gzip_decoded_by_runner": False,
            "graph_accessed": False,
            "source_values_selected": False,
            "ordinary_Z3_semantics_proved": False,
            "physical_radius_proved": False,
            "numeric_error_computed": False,
            "global_closure_proved": False,
            "runtime_model_metadata": None,
            "result": {"z3_row_mapping_capture": _json_safe(source_result)},
            "source_manifest": _json_safe(source),
            "error": None,
        }
        _guard_false_scope(output, "runner success result")
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "manifest_paths_exact": sorted(hashes),
            "leaf_imported": leaf_imported,
            "leaf_callable": call_name,
            "leaf_api_invocations": int(call_name is not None),
            "synthetic_fixture_invocations": 0,
            "source_consumer_invocations": int(call_name is not None),
            "gzip_decoded_by_runner": False,
            "graph_accessed": False,
            "source_values_selected": False,
            "ordinary_Z3_semantics_proved": False,
            "physical_radius_proved": False,
            "numeric_error_computed": False,
            "global_closure_proved": False,
            "runtime_model_metadata": None,
            "result": {"z3_row_mapping_capture": _json_safe(source_result)},
            "source_manifest": _json_safe(source),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }

    result_sha = _write_exclusive(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "source_checkpoint_sha256": source_sha,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "leaf_callable": call_name,
                "leaf_api_invocations": output["leaf_api_invocations"],
                "synthetic_fixture_invocations": 0,
                "source_consumer_invocations": output["source_consumer_invocations"],
                "index": 3,
                "numeric_use_admissible": False,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
