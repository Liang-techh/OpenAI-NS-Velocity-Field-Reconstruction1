"""One frozen raw-authenticated entry for the v2 core P/Dg leaf.

The v1 failure archive is an authenticated raw-only input.  It is intentionally
never decompressed or decoded.  The five successful receipts are decoded only
after the v2 synthetic fixture has passed.  Saved feedback remains Dg input;
the independently integrated compensation receipt remains P provenance only.
The actual v2 leaf is called once and its conditional result is not promoted to
field identity because the live owner is not established by this runner.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import importlib.machinery
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
COMPONENT = "n1_core_p_dg_correlation_v2"
PHASE = COMPONENT + "_pending"
PREVIOUS_COMPONENT = "n1_core_p_dg_correlation_v1"

BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_core_p_dg_correlation_v2.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_core_p_dg_correlation_fixture_v2.py"
RUNNER_REL = "scripts/run_n1_core_p_dg_correlation_v2_once.py"
LEAF_MODULE = "openai_ns_reconstruction.n1_core_p_dg_correlation_v2"
FIXTURE_MODULE = "openai_ns_reconstruction.n1_core_p_dg_correlation_fixture_v2"

RECEIPTS = {
    "KCHI": "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz",
    "FEEDBACK": "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz",
    "COMP": "experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz",
    "LOW": "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz",
    "LOW3": "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz",
}
FAILURE_ARCHIVE = "experiments/root_st073/n1_core_p_dg_correlation_v1_failure.json.gz"
PINNED_RECEIPTS = {
    RECEIPTS["FEEDBACK"]: "0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa",
    RECEIPTS["COMP"]: "c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf",
    RECEIPTS["KCHI"]: "b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0",
    RECEIPTS["LOW3"]: "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6",
    RECEIPTS["LOW"]: "c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9",
}
BASE_INPUTS = frozenset({LEAF_REL, FIXTURE_REL, RUNNER_REL, *RECEIPTS.values(), FAILURE_ARCHIVE})
RAW_ONLY_SUFFIXES = ("_failure.json.gz", "_failure.jsonl.gz")
FORBIDDEN_PATH_PARTS = (
    "docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json",
    "Rm_loop_density_identity.json",
)
FORBIDDEN_MODULE_BASES = frozenset(
    {
        "n1_core_p_dg_correlation_v1",
        "n1_core_p_dg_correlation_fixture_v1",
        "n1_core_kchi_feedback_v2",
        "n1_core_kchi_compensation_v2",
        "n1_kchi_physical_compensation_v2",
        "n1_core_kchi_first_v1",
        "n1_core_low_word_v3",
        "n1_core_word1_z2_v2",
        "n1_local_p_dg_correlated_pair_v1",
        "n1_local_p_dg_correlated_pair_fixture_v1",
    }
)
FORBIDDEN_FILE_NAMES = frozenset(name + ".py" for name in FORBIDDEN_MODULE_BASES)
REQUIRED_ASSERTIONS = (
    "same_actual_g5",
    "same_chi_Z",
    "same_green_K0",
    "same_complete_kernel",
    "same_physical_m15_measure",
    "same_leading_F",
    "low3_W0_same_K0_g5",
    "kchi_Dg_identity",
    "ordinary_Z1_chi_Z",
    "ordinary_Z1_g5_Z",
    "same_geometry",
    "same_labels",
    "same_parent_hashes",
    "exact_source_id_paths",
    "exact_cell_paths",
    "exact_outputs",
)


class RunnerInputError(ValueError):
    """Raised when the frozen v2 one-entry contract is not satisfied."""


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
    raw_interval = getattr(value, "_mpi_", None)
    if raw_interval is not None:
        return {
            "lower_exact_mpf_tuple": list(raw_interval[0]),
            "upper_exact_mpf_tuple": list(raw_interval[1]),
            "directed": True,
        }
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
            return sorted((_json_safe(item, active) for item in value), key=repr)
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


def _decode_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        return _decode_json(path.read_bytes(), label)
    except OSError as exc:
        raise RunnerInputError(label + " is unreadable") from exc


def _canonical_repo_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    _require(not any(part in normalized for part in FORBIDDEN_PATH_PARTS), label + " is forbidden")
    relative = Path(normalized)
    _require(not relative.is_absolute() and ".." not in relative.parts, label + " must be repository-relative")
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(), label + " is missing or escapes repository")
    return key, path


def _brief_file_path(raw_path: Any, label: str) -> Path:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    _require(not any(part in normalized for part in FORBIDDEN_PATH_PARTS), label + " is forbidden")
    candidate = Path(normalized)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    candidate = candidate.resolve()
    _require(candidate.is_file(), label + " is missing")
    return candidate


def _sha_text(value: Any, label: str) -> str:
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


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


def _brief_checkpoint(brief: Mapping[str, Any], record: Mapping[str, Any]) -> Mapping[str, Any]:
    for candidate in (
        record.get("source_checkpoint"),
        brief.get("source_checkpoint"),
        brief.get(COMPONENT + "_source_checkpoint"),
    ):
        if isinstance(candidate, Mapping):
            return candidate
    raise RunnerInputError("source checkpoint is missing")


def _manifest_hashes(source: Mapping[str, Any]) -> tuple[dict[str, str], list[str]]:
    files = source.get("files")
    _require(isinstance(files, list) and files, "source manifest files are missing")
    hashes: dict[str, str] = {}
    raw_only: list[str] = []
    for index, row in enumerate(files):
        _require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        key, path = _canonical_repo_path(row.get("path"), f"source.files[{index}]")
        _require(key not in hashes, "duplicate source path " + key)
        expected = _sha_text(row.get("sha256"), key)
        actual = _sha_file(path)
        _require(actual == expected, "raw SHA mismatch: " + key)
        if key in PINNED_RECEIPTS:
            _require(actual == PINNED_RECEIPTS[key], "pinned receipt changed: " + key)
        if key.endswith(".json.gz"):
            _require(
                row.get("raw_sha256_before_decode", row.get("authenticate_raw_before_decode")) is True,
                key + " is not marked raw-authenticated-before-decode",
            )
        if key.endswith(RAW_ONLY_SUFFIXES) or key == FAILURE_ARCHIVE:
            raw_only.append(key)
        hashes[key] = actual
    _require(BASE_INPUTS.issubset(hashes), "manifest is missing a required input")
    _require(FAILURE_ARCHIVE in raw_only, "v1 failure archive is not raw-only")
    return hashes, raw_only


def _preflight() -> tuple[dict[str, Any], dict[str, str], list[str], str, str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(BRIEF.is_file(), "resume brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = _decode_json(brief_raw, "resume brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending brief one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending brief reentry guard missing")
    if "goal_active" in record:
        _require(record.get("goal_active") is True, "pending brief goal is not active")
    if "full_goal_complete" in record:
        _require(record.get("full_goal_complete") is False, "pending brief is already complete")
    _require(record.get("previous_component", PREVIOUS_COMPONENT) == PREVIOUS_COMPONENT, "pending brief previous component mismatch")
    checkpoint = _brief_checkpoint(brief, record)
    checkpoint_path = _brief_file_path(checkpoint.get("path"), "brief source checkpoint")
    _require(checkpoint_path == SOURCE.resolve(), "brief points to another source")
    source_sha = _sha_file(SOURCE)
    _require(_sha_text(checkpoint.get("sha256"), "source checkpoint") == source_sha, "source checkpoint SHA mismatch")
    preimage_sha = ""
    if checkpoint.get("preimage_path"):
        preimage = _brief_file_path(checkpoint.get("preimage_path"), "brief preimage")
        preimage_sha = _sha_file(preimage)
        _require(_sha_text(checkpoint.get("preimage_sha256"), "brief preimage") == preimage_sha, "brief preimage SHA mismatch")
    source = _load_json(SOURCE, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")
    _require(source.get("previous_component") == PREVIOUS_COMPONENT, "source manifest previous component mismatch")
    contract = source.get("input_contract")
    if isinstance(contract, Mapping):
        _require("build(c, dg_physical, certificate, source_receipts)" in str(contract.get("leaf_api", "")), "leaf API contract mismatch")
        _require("run_fixture" in str(contract.get("fixture_api", "")), "fixture API contract mismatch")
    hashes, raw_only = _manifest_hashes(source)
    return source, hashes, raw_only, source_sha, _sha_bytes(brief_raw), preimage_sha


def _forbidden_module_names() -> list[str]:
    found: list[str] = []
    for name, module in sys.modules.items():
        base = name.rsplit(".", 1)[-1]
        module_file = getattr(module, "__file__", None)
        file_name = Path(module_file).name if isinstance(module_file, str) else ""
        if base in FORBIDDEN_MODULE_BASES or file_name in FORBIDDEN_FILE_NAMES:
            found.append(name)
    return sorted(set(found))


def _load_new_modules() -> tuple[Any, Any]:
    _require(not _forbidden_module_names(), "historical production module already loaded")
    package_name = "openai_ns_reconstruction"
    package_dir = SRC / package_name
    previous_package = sys.modules.get(package_name)
    previous_leaf = sys.modules.get(LEAF_MODULE)
    previous_fixture = sys.modules.get(FIXTURE_MODULE)
    namespace = ModuleType(package_name)
    namespace.__path__ = [str(package_dir)]
    namespace.__package__ = package_name
    namespace.__spec__ = importlib.machinery.ModuleSpec(package_name, loader=None, is_package=True)
    sys.modules[package_name] = namespace
    old_path = list(sys.path)
    source_root = str(SRC)
    if source_root not in sys.path:
        sys.path.insert(0, source_root)
    try:
        leaf_spec = importlib.util.spec_from_file_location(LEAF_MODULE, ROOT / LEAF_REL)
        fixture_spec = importlib.util.spec_from_file_location(FIXTURE_MODULE, ROOT / FIXTURE_REL)
        _require(leaf_spec is not None and leaf_spec.loader is not None, "v2 leaf module spec is unavailable")
        _require(fixture_spec is not None and fixture_spec.loader is not None, "v2 fixture module spec is unavailable")
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
    _require(not _forbidden_module_names(), "v2 leaf/fixture imported a historical producer")
    return leaf, fixture


def _decode_intervals(context: MPIntervalContext, value: Any) -> Any:
    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            _require(has_lower and has_upper, "incomplete directed interval endpoints")
            lower = tuple(value["lower_exact_mpf_tuple"])
            upper = tuple(value["upper_exact_mpf_tuple"])
            _require(len(lower) == len(upper) == 4 and all(type(item) is int for item in lower + upper), "invalid exact directed interval endpoints")
            result = context.make_mpf((lower, upper))
            _require(result.a <= result.b, "reversed directed interval endpoints")
            return result
        return {key: _decode_intervals(context, item) for key, item in value.items()}
    if isinstance(value, list):
        return [_decode_intervals(context, item) for item in value]
    if isinstance(value, tuple):
        return [_decode_intervals(context, item) for item in value]
    return value


def _decode_receipt(raw: Mapping[str, bytes], hashes: Mapping[str, str], label: str, context: MPIntervalContext) -> dict[str, Any]:
    path = RECEIPTS[label]
    try:
        value = json.loads(gzip.decompress(raw[path]))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError("saved receipt is not valid gzip JSON: " + path) from exc
    _require(isinstance(value, Mapping), "saved receipt is not an object: " + path)
    result = value.get("result")
    _require(isinstance(result, Mapping) and result.get("status") == "returned", "saved receipt is not returned: " + path)
    source_manifest = value.get("source")
    _require(isinstance(source_manifest, Mapping), "saved receipt source manifest is missing: " + path)
    output = result.get("output")
    _require(output is not None, "saved receipt output is missing: " + path)
    return {
        "path": path,
        "raw_sha256": hashes[path],
        "raw_hash_authenticated": True,
        "decoded": True,
        "status": "returned",
        "source_manifest": dict(source_manifest),
        "decoded_output": _decode_intervals(context, output),
        "raw_receipt": value,
    }


def _path_get(value: Any, *keys: Any) -> Any:
    current = value
    for key in keys:
        if not isinstance(current, Mapping):
            return None
        if key in current:
            current = current[key]
        elif str(key) in current:
            current = current[str(key)]
        else:
            return None
    return current


def _flatten_saved_row(value: Any) -> list[Any]:
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, Mapping):
        flattened: list[Any] = []
        for item in value.values():
            if isinstance(item, (list, tuple)):
                flattened.extend(item)
        return flattened
    return []


def _find_saved_row(output: Mapping[str, Any], candidates: Sequence[tuple[Any, ...]], label: str) -> tuple[tuple[Any, ...], list[Any]]:
    for path in candidates:
        value = _path_get(output, *path)
        flattened = _flatten_saved_row(value)
        if flattened:
            return path, flattened
    raise RunnerInputError(label)


def _feedback_rows(output: Mapping[str, Any], order: int) -> tuple[tuple[Any, ...], list[Any], tuple[Any, ...], list[Any]]:
    signed_candidates = (
        ("m15", "swirl", str(order), "W0_chi"),
        ("m15", "swirl", order, "W0_chi"),
        ("m15", "swirl_by_order", str(order), "W0_chi"),
        ("m15", "swirl_by_order", order, "W0_chi"),
        ("m15", "swirl_by_branch", "W0_chi", str(order)),
        ("m15", "swirl_by_branch", "W0_chi", order),
        ("m15", "branches", "W0_chi", str(order), "swirl"),
        ("m15", "branches", "W0_chi", order, "swirl"),
    )
    error_candidates = (
        ("m15", "swirl_error", str(order), "W0_chi"),
        ("m15", "swirl_error", order, "W0_chi"),
        ("m15", "swirl_errors_by_order", str(order), "W0_chi"),
        ("m15", "swirl_errors_by_order", order, "W0_chi"),
        ("m15", "symmetric_error_swirl", str(order), "W0_chi"),
        ("m15", "symmetric_error_swirl", order, "W0_chi"),
        ("m15", "error_sectors", "W0_chi", str(order), "swirl"),
        ("m15", "error_sectors", "W0_chi", order, "swirl"),
        ("m15", "error_sectors", "W0_chi", str(order), "swirl_error"),
        ("m15", "error_sectors", "W0_chi", order, "swirl_error"),
        ("m15", "branches", "W0_chi", str(order), "swirl_error"),
        ("m15", "branches", "W0_chi", order, "swirl_error"),
    )
    signed_path, signed = _find_saved_row(output, signed_candidates, f"feedback saved W0_chi swirl ordinary Z{order} row is missing")
    error_path, errors = _find_saved_row(output, error_candidates, f"feedback saved W0_chi swirl_error ordinary Z{order} row is missing")
    return signed_path, signed, error_path, errors


def _source_items(term: Mapping[str, Any], origin: Mapping[str, Any], plural: str, singular: str) -> list[Any]:
    values: list[Any] = []
    for container in (term, origin):
        for key in (plural, singular):
            value = container.get(key)
            if isinstance(value, (list, tuple, set, frozenset)):
                values.extend(value)
            elif value is not None:
                values.append(value)
    unique: list[Any] = []
    for value in values:
        if not any(value == prior for prior in unique):
            unique.append(value)
    return unique


def _dg_packet(feedback_output: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    rows: dict[str, list[dict[str, Any]]] = {}
    index: dict[str, Any] = {}
    row_paths: dict[str, Any] = {}
    for order in (0, 1):
        signed_path, signed_terms, error_path, error_terms = _feedback_rows(feedback_output, order)
        z_label = f"Z{order}"
        row_paths[z_label] = {
            "swirl": ".".join(str(item) for item in signed_path),
            "swirl_error": ".".join(str(item) for item in error_path),
        }
        converted: list[dict[str, Any]] = []
        saved_routes = ((False, signed_path, signed_terms), (True, error_path, error_terms))
        for error_route, path, saved_terms in saved_routes:
            for index_in_row, saved in enumerate(saved_terms):
                _require(isinstance(saved, Mapping), f"feedback {z_label} Dg term is not an object")
                term = copy.deepcopy(dict(saved))
                origin = term.get("origin") if isinstance(term.get("origin"), Mapping) else {}
                source_ids = _source_items(term, origin, "source_ids", "source_id")
                source_paths = _source_items(term, origin, "source_paths", "source_path")
                if not source_ids and source_paths:
                    source_ids = list(source_paths)
                    term["source_id_derived_from_saved_path"] = True
                _require(source_ids, f"feedback {z_label} Dg term has no saved source reference")
                if not source_paths:
                    source_paths = [f"{RECEIPTS['FEEDBACK']}:result.output.{'.'.join(str(item) for item in path)}.{index_in_row}"]
                parent_id = term.get("parent_id")
                if parent_id is None:
                    parent_id = f"saved-Dg:{z_label}:{'error' if error_route else 'swirl'}:{index_in_row}"
                if error_route:
                    term["error_only"] = True
                    term["symmetric_error"] = True
                    term["uncertainty"] = True
                else:
                    if "error_only" not in term and "error_only" in origin:
                        term["error_only"] = origin["error_only"]
                    if "symmetric_error" not in term and "symmetric_error" in origin:
                        term["symmetric_error"] = origin["symmetric_error"]
                term.update(
                    {
                        "parent_id": parent_id,
                        "source_ids": source_ids,
                        "source_paths": source_paths,
                        "formal_powers": list(term.get("formal_powers", term.get("powers", []))),
                        "saved_receipt_path": RECEIPTS["FEEDBACK"],
                        "saved_output_path": ["result", "output", *path, index_in_row],
                    }
                )
                converted.append(term)
                index[str(parent_id)] = {
                    "receipt": RECEIPTS["FEEDBACK"],
                    "path": ["result", "output", *path, index_in_row],
                    "source_ids": source_ids,
                    "source_paths": source_paths,
                    "error_only": bool(term.get("error_only", False)),
                }
        rows[z_label] = converted
    return {"physical_level": "m15", "J1": rows}, {"row_paths": row_paths, "term_index": index}


def _owner_records(packet: Mapping[str, Any]) -> dict[str, Any]:
    owners: dict[str, Any] = {}
    for z_label in ("Z0", "Z1"):
        for index, term in enumerate(packet["J1"][z_label]):
            parent = str(term["parent_id"])
            owners[parent] = {
                "parent_id": term["parent_id"],
                "owner": "saved Dg physical",
                "passed": False,
                "proved": False,
                "source_ids": list(term["source_ids"]),
                "source_paths": list(term["source_paths"]),
                "parent_hash": "saved receipt only; live owner not proved",
                "geometry": "saved feedback m15 row only",
                "labels": {"ordinary_Z": z_label, "row_index": index},
                "exact_cell_path": ".".join(str(item) for item in term["saved_output_path"]),
                "output_path": ".".join(str(item) for item in term["saved_output_path"]),
            }
    return owners


def _saved_exact_assertion(outputs: Sequence[Mapping[str, Any]], names: Sequence[str]) -> bool:
    wanted = set(names)
    return any(output.get(name) is True for output in outputs for name in wanted)


def _build_certificate(receipts: Mapping[str, Mapping[str, Any]], packet: Mapping[str, Any], dg_index: Mapping[str, Any]) -> dict[str, Any]:
    outputs = [entry["decoded_output"] for entry in receipts.values()]
    certificate: dict[str, Any] = {
        "source_receipts": {
            label: {key: value for key, value in entry.items() if key != "raw_receipt"}
            for label, entry in receipts.items()
        },
        "owner_by_parent": _owner_records(packet),
        "Dg_source_refs": {"receipt": RECEIPTS["FEEDBACK"], "output": "result.output.m15.branches.W0_chi", "row_paths": dg_index["row_paths"], "term_index": dg_index["term_index"]},
        "P_source_refs": {"receipt": RECEIPTS["COMP"], "output": "result.output.potential_compensation_m15", "independent": True, "not_substituted_for_Dg": True},
        "P_sign": 2,
        "certificate_scope": "saved receipt references only; live same-owner proof is not asserted",
    }
    exact_specs = {
        "same_actual_g5": ("same_actual_g5", "same_g5_function", "same_g5_source_function", "g5_exact_owner_verified"),
        "same_chi_Z": ("same_chi_Z", "same_chi_function", "same_chi_source_function", "chi_Z_same_source"),
        "same_green_K0": ("same_green_K0", "same_K0_function", "K0_exact_same_source", "same_zero_inlet_green"),
        "same_complete_kernel": ("same_complete_0_2_kernel", "same_complete_[0,2]_kernel", "complete_kernel_0_2_same", "same_complete_kernel"),
        "same_physical_m15_measure": ("same_physical_m15_measure", "same_m15_measure", "physical_m15_measure_same", "same_physical_measure_cutoff"),
        "same_leading_F": ("same_leading_F", "same_leading_F0", "leading_F_same_source", "same_leading_field"),
        "low3_W0_same_K0_g5": ("low3_W0_same_K0_g5", "low3_W0_exact_K0_g5", "W0_same_K0_g5", "low3_exact_function_owner"),
        "kchi_Dg_identity": ("kchi_Dg_minus2_chi_K0_K0g5", "kchi_Dg_same_minus2_chi_K0_K0g5", "Dg_same_minus2_chi_K0_K0g5", "kchi_Dg_function_identity"),
        "ordinary_Z1_chi_Z": ("ordinary_Z1_chi_Z_retained", "ordinary_Z1_retains_chi_Z", "Z1_chi_Z_present"),
        "ordinary_Z1_g5_Z": ("ordinary_Z1_g5_Z_retained", "ordinary_Z1_retains_g5_Z", "Z1_g5_Z_present"),
        "same_geometry": ("same_geometry", "geometry_same_source", "physical_geometry_same"),
        "same_labels": ("same_labels", "source_labels_same", "physical_labels_same"),
        "same_parent_hashes": ("same_parent_hashes", "parent_hashes_same", "all_parent_hashes_match"),
        "exact_source_id_paths": ("exact_source_id_paths", "source_id_paths_verified", "source_paths_exact"),
        "exact_cell_paths": ("exact_cell_paths", "cell_paths_verified", "source_cell_paths_exact"),
        "exact_outputs": ("exact_outputs", "outputs_verified", "source_outputs_exact"),
    }
    for key, names in exact_specs.items():
        certificate[key] = _saved_exact_assertion(outputs, names)
    certificate["saved_exact_assertions"] = {key: certificate[key] for key in exact_specs}
    return certificate


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not an object")
    result = dict(value)
    _require(result.get("schema") == "n1_core_p_dg_correlation_fixture_v2", "v2 fixture schema mismatch")
    _require(result.get("passed") is True, "v2 synthetic fixture did not pass")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()), "v2 fixture checks did not pass")
    _require(result.get("synthetic_only") is True, "v2 fixture is not synthetic-only")
    _require(result.get("implementation_imported") is False, "v2 fixture imported implementation directly")
    _require(result.get("old_code_called") is False, "v2 fixture called old code")
    return result


def _summarise_compensation(value: Mapping[str, Any]) -> dict[str, Any]:
    branches = value.get("potential_compensation_m15")
    return {
        "schema": value.get("schema"),
        "branch_names": sorted(str(key) for key in branches) if isinstance(branches, Mapping) else [],
        "ordinary_Z_orders": value.get("ordinary_Z_orders"),
        "formula": value.get("P_formula"),
        "independent_data_retained": True,
        "used_as_Dg": False,
        "used_as_per_sector_P_alias": False,
    }


def _consume(hashes: Mapping[str, str], raw_only: Sequence[str], source_sha: str, brief_sha: str, preimage_sha: str, progress: dict[str, Any]) -> dict[str, Any]:
    raw: dict[str, bytes] = {}
    for path, digest in hashes.items():
        data = (ROOT / path).read_bytes()
        _require(_sha_bytes(data) == digest, "input changed after entry: " + path)
        raw[path] = data
    progress["raw_sha_reauthenticated"] = True
    progress["failure_archives"] = {
        path: {"raw_sha256": hashes[path], "raw_hash_authenticated": True, "decoded": False, "decompressed": False}
        for path in raw_only
    }

    leaf, fixture = _load_new_modules()
    progress["historical_modules_imported"] = _forbidden_module_names()
    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "v2 fixture must export run_fixture")
    context = MPIntervalContext()
    context.prec = 1800
    progress["fixture_calls"] = 1
    fixture_result = _validate_fixture(run_fixture(leaf.build, context))
    progress["fixture"] = _json_safe(fixture_result)
    _require(not _forbidden_module_names(), "v2 fixture imported a historical producer")

    decoded = {label: _decode_receipt(raw, hashes, label, context) for label in ("FEEDBACK", "COMP", "KCHI", "LOW3", "LOW")}
    progress["decoded_success_receipts"] = list(decoded)
    feedback_output = decoded["FEEDBACK"]["decoded_output"]
    _require(isinstance(feedback_output, Mapping), "feedback output is not an object")
    dg_packet, dg_index = _dg_packet(feedback_output)
    evidence = {label: {key: value for key, value in entry.items() if key != "raw_receipt"} for label, entry in decoded.items()}
    certificate = _build_certificate(decoded, dg_packet, dg_index)
    progress["certificate"] = _json_safe(certificate)
    compensation_output = decoded["COMP"]["decoded_output"]
    _require(isinstance(compensation_output, Mapping), "compensation output is not an object")
    progress["independent_compensation"] = _summarise_compensation(compensation_output)
    progress["source_provenance"] = {
        "same_family_is_not_function_identity": True,
        "same_owner_proved": False,
        "same_K0_proved": False,
        "same_m15_measure_proved": False,
        "independent_P_not_substituted_for_Dg": True,
    }

    build = getattr(leaf, "build", None)
    _require(callable(build), "v2 leaf must export build")
    progress["actual_leaf_calls"] = 1
    actual = build(context, dg_packet, certificate, evidence)
    encoded_actual = _json_safe(actual)
    progress["actual_result"] = encoded_actual
    _require(isinstance(actual, Mapping), "v2 leaf result is not an object")
    _require(not _forbidden_module_names(), "v2 leaf imported a historical producer")

    _require(actual.get("status") in {"conditional", "conditional_on_caller_assertions"}, "v2 leaf did not remain conditional")
    _require(
        actual.get("passed") is False
        and actual.get("correlation_proved") is False
        and actual.get("actual_source_verified") is False
        and actual.get("identity_verified") is False,
        "v2 conditional result claimed field identity",
    )
    if actual.get("status") == "conditional_on_caller_assertions":
        _require(isinstance(actual.get("functional_pairs"), list) and len(actual["functional_pairs"]) == 1, "v2 conditional whole-function pair is missing")
    else:
        _require(not actual.get("functional_pairs"), "v2 failed-evidence branch manufactured a functional pair")
    p_rows = actual.get("P", {})
    _require(isinstance(p_rows, Mapping) and not p_rows.get("Z0") and not p_rows.get("Z1"), "v2 result manufactured per-sector P rows")
    dg_rows = actual.get("Dg", {})
    _require(isinstance(dg_rows, Mapping) and dg_rows.get("Z0") is not None and dg_rows.get("Z1") is not None, "v2 result dropped saved Dg rows")
    return {
        "schema": COMPONENT + "_runner_result_v2",
        "status": "returned",
        "passed": False,
        "leaf_status": actual.get("status"),
        "leaf_correlation_proved": actual.get("correlation_proved") is True,
        "leaf_actual_source_verified": actual.get("actual_source_verified") is True,
        "conditional_not_field_identity": True,
        "fixture": fixture_result,
        "actual_result": encoded_actual,
        "source_receipts": {
            label: {"path": entry["path"], "raw_sha256": entry["raw_sha256"], "raw_hash_authenticated": entry["raw_hash_authenticated"], "decoded": entry["decoded"], "status": entry["status"], "source_manifest": entry["source_manifest"]}
            for label, entry in decoded.items()
        },
        "Dg_source_index": dg_index,
        "independent_compensation": progress["independent_compensation"],
        "source_provenance": progress["source_provenance"],
        "failure_archives": progress["failure_archives"],
        "input_hashes": dict(hashes),
        "source_sha256": source_sha,
        "brief_sha256": brief_sha,
        "preimage_sha256": preimage_sha,
        "actual_leaf_calls": 1,
        "ancestor_math_calls": 0,
        "historical_producer_imported": False,
        "old_errors_removed": False,
        "interval_cancellation_performed": False,
        "independent_P_substituted_for_Dg": False,
        "ledger_written": False,
        "terminal_closure_proved": False,
    }


def main() -> int:
    source, hashes, raw_only, source_sha, brief_sha, preimage_sha = _preflight()
    del source
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": COMPONENT + "_entered_v2",
            "component": COMPONENT,
            "phase": PHASE,
            "previous_component": PREVIOUS_COMPONENT,
            "entered": True,
            "source_manifest_sha256": source_sha,
            "brief_sha256": brief_sha,
            "preimage_sha256": preimage_sha,
            "input_hashes": hashes,
            "raw_sha_authenticated_before_entry": True,
            "receipt_decode_started": False,
            "failure_archives_raw_only": {path: True for path in raw_only},
            "actual_leaf_calls": 0,
            "ancestor_math_calls": 0,
            "historical_producer_imported": False,
        },
    )
    progress: dict[str, Any] = {"actual_leaf_calls": 0, "fixture_calls": 0, "raw_sha_reauthenticated": False, "decoded_success_receipts": []}
    try:
        nested = _consume(hashes, raw_only, source_sha, brief_sha, preimage_sha, progress)
        output = {"component": COMPONENT, "entry_sha256": entry_sha, **nested, "error": None, "traceback": None}
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_runner_result_v2",
            "component": COMPONENT,
            "status": "failed",
            "passed": False,
            "entry_sha256": entry_sha,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "preimage_sha256": preimage_sha,
            "input_hashes": dict(hashes),
            "progress": _json_safe(progress),
            "failure_archives": progress.get("failure_archives", {path: {"raw_sha256": hashes[path], "raw_hash_authenticated": True, "decoded": False, "decompressed": False} for path in raw_only}),
            "actual_leaf_calls": int(progress.get("actual_leaf_calls", 0)),
            "ancestor_math_calls": 0,
            "historical_producer_imported": bool(_forbidden_module_names()),
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        }
    result_sha = _write_exclusive(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"], "result_sha256": result_sha}, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
