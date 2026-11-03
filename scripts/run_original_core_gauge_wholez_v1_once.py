"""Run the new original-core gauge whole-Z leaf exactly once.

The runner owns the one-entry boundary for ``original_core_gauge_wholez_v1``.
The root agent freezes ``tmp/original_core_gauge_wholez_v1/source.json`` and
the global continuation brief before this runner is invoked.  The manifest
contains this runner, the new leaf, its new fixture, the relevant original
source texts, the prior source-review documents, and the required provenance
records.  The global continuation brief and its pending record are
authenticated separately rather than listed in ``source.json.files``; this
avoids a manifest-SHA/pending-brief self-reference.  The runner writes an
exclusive entry receipt, runs the independent fixture once, and then uses one
new leaf owner invocation to evaluate four contiguous cells covering
``[-1, 1]``.

No archived gzip member is decoded, no old producer or old mathematical
fixture is imported, and no old integral or four-cell result is recomputed.
The result is deliberately not an NS closure, a live-owner certificate, or
an equality claim obtained by treating the upper delta bound as a value.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import inspect
import json
import os
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "original_core_gauge_wholez_v1"
PHASE = COMPONENT + "_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"
GLOBAL_BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
GLOBAL_BRIEF_REL = "../../RESUME_BRIEF_2026_10_16.json"

LEAF = "src/openai_ns_reconstruction/original_core_gauge_wholez_v1.py"
FIXTURE = "src/openai_ns_reconstruction/original_core_gauge_wholez_v1_fixture.py"
RUNNER = "scripts/run_original_core_gauge_wholez_v1_once.py"

# These are source/check files, not producer inputs.  The runner hashes them
# for provenance only; it never imports them and never executes their checks.
ORIGINAL_SOURCE_FILES = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source_check.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_core_hierarchy_source.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_core_hierarchy_source_check.py",
)

PRIOR_DOCUMENTS = (
    "docs/CURRENT_N1_CORE_GAUGE_EXIT_2026_10_30.md",
    "docs/CURRENT_ORIGINAL_N1_CORE_FIXED_POINT_FUNCTIONS_2026_10_13.md",
    "docs/CURRENT_ORIGINAL_N1_ANALYTIC_CORE_SOURCE_2026_10_11.md",
    "docs/CURRENT_N1_RAW_SHEAR_CHI_INTERVAL_2026_11_02.md",
    "docs/CURRENT_N1_COMMON_SOURCE_BRIDGE_2026_11_02.md",
)

PROVENANCE_FILES = (
    "docs/CURRENT_ORIGINAL_CORE_LOG_AMPLITUDE_2026_10_15.md",
    "experiments/root_st073/lei_ren_part1_paper_compliant_core_transfer.json",
)

STATIC_FROZEN_INPUTS = (
    (LEAF, FIXTURE, RUNNER)
    + ORIGINAL_SOURCE_FILES
    + PRIOR_DOCUMENTS
    + PROVENANCE_FILES
)

J = Fraction(1, 5721363714816000000000)
DELTA_INTERVAL = (Fraction(0), Fraction(1, 10**200))
ANCHOR_DOMAIN = (-J, Fraction(0))
CELLS = (
    (Fraction(-1), Fraction(-1, 2)),
    (Fraction(-1, 2), Fraction(0)),
    (Fraction(0), Fraction(1, 2)),
    (Fraction(1, 2), Fraction(1)),
)

LEAF_API = "ActualCoreGaugeWholeZ(source=GaugeSource(...)).evaluate_cell(left, right) for four cells"
FIXTURE_API = "run_fixture(leaf)"
SCHEMA = COMPONENT + "_result_v1"
MANIFEST_SCHEMA = COMPONENT + "_manifest_v1"

# New code must be self-contained for this one-entry experiment.  In
# particular it must not import any of the historical original_* producers.
FORBIDDEN_IMPORT_PREFIXES = (
    "openai_ns_reconstruction.original_",
    "experiments.",
    "scripts.",
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


def require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_bytes(value: Any) -> bytes:
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


def json_safe(value: Any, seen: set[int] | None = None) -> Any:
    """Serialize intervals and exact values without selecting a midpoint."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Fraction):
        return {
            "kind": "exact_fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
            "source_value_used": False,
        }
    if isinstance(value, float):
        require(value == value and abs(value) != float("inf"), "non-finite float in result")
        return value
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"kind": "bytes", "length": len(value), "sha256": sha_bytes(value)}
    if hasattr(value, "_mpi_"):
        try:
            lower, upper = value._mpi_
            return {
                "kind": "interval_enclosure",
                "lower_exact_mpf_tuple": list(tuple(lower)),
                "upper_exact_mpf_tuple": list(tuple(upper)),
                "source_value_used": False,
                "midpoint_used_as_function_value": False,
                "delta_upper_used_as_value": False,
            }
        except Exception:
            return {"kind": "interval_enclosure", "repr": repr(value), "source_value_used": False}
    if hasattr(value, "_mpf_"):
        return {
            "kind": "exact_mpf",
            "exact_mpf_tuple": list(tuple(value._mpf_)),
            "source_value_used": False,
        }
    if seen is None:
        seen = set()
    identity = id(value)
    if identity in seen:
        return {"kind": "cycle"}
    seen.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): json_safe(item, seen) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(item, seen) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((json_safe(item, seen) for item in value), key=repr)
        # A leaf result must be serializable, but preserving an opaque repr is
        # safer than losing the complete traceback when a new packet is added.
        return {"kind": "opaque", "type": type(value).__name__, "repr": repr(value)}
    finally:
        seen.discard(identity)


def save_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json_bytes(json_safe(value))
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
    return sha_bytes(raw)


def repository_path(relative: str) -> Path:
    path = Path(relative.replace("/", os.sep))
    require(not path.is_absolute() and ".." not in path.parts, "non-relative input path: " + relative)
    resolved = (ROOT / path).resolve()
    require(resolved.is_relative_to(ROOT.resolve()), "input escapes repository: " + relative)
    require(resolved.is_file(), "missing frozen input: " + relative)
    return resolved


def load_json_bytes(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("utf-8"),
            parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite JSON") from exc
    require(isinstance(value, Mapping), label + " must be a JSON object")
    return dict(value)


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        return load_json_bytes(path.read_bytes(), label)
    except OSError as exc:
        raise RunnerInputError(label + " cannot be read") from exc


def resolve_pending_path(raw_path: Any) -> Path:
    require(isinstance(raw_path, str) and raw_path, "global pending brief path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        from_global = (GLOBAL_BRIEF.parent / candidate).resolve()
        from_repo = (ROOT / candidate).resolve()
        candidate = from_global if from_global.is_file() else from_repo
    else:
        candidate = candidate.resolve()
    require(candidate.is_file(), "global pending brief file missing: " + str(candidate))
    return candidate


def canonical_key(path: Path) -> str:
    path = path.resolve()
    if path == GLOBAL_BRIEF.resolve():
        return GLOBAL_BRIEF_REL
    try:
        return path.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def frozen_path(key: str, pending_path: Path) -> Path:
    if key == GLOBAL_BRIEF_REL:
        return GLOBAL_BRIEF.resolve()
    candidate = Path(key)
    if candidate.is_absolute():
        resolved = candidate.resolve()
        if resolved == GLOBAL_BRIEF.resolve():
            return resolved
        require(resolved == pending_path.resolve(), "external frozen path is not the pending brief")
        require(resolved.is_file(), "external frozen input missing: " + key)
        return resolved
    return repository_path(key)


def read_global_pending() -> tuple[dict[str, Any], bytes, Path, dict[str, Any], bytes]:
    require(GLOBAL_BRIEF.is_file(), "global continuation brief is missing")
    global_raw = GLOBAL_BRIEF.read_bytes()
    global_doc = load_json_bytes(global_raw, "global continuation brief")
    pointer = global_doc.get("pending_brief")
    require(isinstance(pointer, Mapping), "global continuation brief pending pointer missing")
    pending_path = resolve_pending_path(pointer.get("path"))
    pending_raw = pending_path.read_bytes()
    expected = pointer.get("sha256")
    require(isinstance(expected, str) and sha_bytes(pending_raw) == expected.lower(), "global pending brief SHA mismatch")
    pending = load_json_bytes(pending_raw, "global pending brief")
    return global_doc, global_raw, pending_path, pending, pending_raw


def preflight_frozen() -> tuple[dict[str, Any], dict[str, str], str, str, str, Path]:
    """Read root-frozen manifest and global pending brief; never create either."""

    require(SOURCE.is_file(), "root-frozen source manifest is missing")
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    source_raw = SOURCE.read_bytes()
    source = load_json_bytes(source_raw, "source manifest")
    source_sha = sha_bytes(source_raw)
    global_doc, global_raw, pending_path, pending, pending_raw = read_global_pending()
    global_sha = sha_bytes(global_raw)
    pending_sha = sha_bytes(pending_raw)
    require(source.get("schema") == MANIFEST_SCHEMA, "source manifest schema mismatch")
    require(source.get("component") == COMPONENT and source.get("phase") == PHASE, "source component mismatch")
    require(source.get("one_entry_only") is True and source.get("do_not_reenter") is True, "source entry guard missing")
    require(source.get("source_manifest_self_hash_forbidden") is True, "source manifest self-hash guard missing")
    require(source.get("root_frozen_before_runner") is True, "root-freeze guard missing")
    require(pending.get("component") == COMPONENT and pending.get("phase") == PHASE, "global pending component mismatch")
    require(pending.get("goal_active") is True and pending.get("full_goal_complete") is False, "global pending goal guard failed")
    require(pending.get("one_entry_only") is True and pending.get("do_not_reenter") is True, "global pending entry guard failed")
    checkpoint = pending.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "global pending source checkpoint missing")
    checkpoint_path = Path(str(checkpoint.get("path")))
    if not checkpoint_path.is_absolute():
        checkpoint_path = (ROOT / checkpoint_path).resolve()
    require(checkpoint_path == SOURCE.resolve(), "global pending source checkpoint path mismatch")
    require(checkpoint.get("sha256") == source_sha, "global pending source checkpoint hash mismatch")
    contract = source.get("input_contract")
    require(isinstance(contract, Mapping), "input contract missing")
    require(contract.get("leaf_api") == LEAF_API and contract.get("fixture_api") == FIXTURE_API, "API contract mismatch")
    require(contract.get("whole_Z_cell_count") == 4 and contract.get("ordinary_Z_orders") == [0, 1], "whole-Z contract mismatch")
    policy = source.get("source_policy")
    require(isinstance(policy, Mapping), "source policy missing")
    for flag in ("authenticate_all_frozen_bytes_before_import", "briefs_separately_authenticated", "no_gzip_decode", "no_old_producer_import", "no_old_fixture_import", "fixture_once_before_leaf", "leaf_once_for_four_cells", "exclusive_entry_and_result", "save_full_traceback_on_failure"):
        require(policy.get(flag) is True, "source policy guard missing: " + flag)
    rows = source.get("files")
    require(isinstance(rows, list), "source manifest files missing")
    expected_paths = {repository_path(path).resolve() for path in STATIC_FROZEN_INPUTS}
    separately_authenticated_briefs = {GLOBAL_BRIEF.resolve(), pending_path.resolve()}
    hashes: dict[str, str] = {}
    resolved_rows: set[Path] = set()
    for row in rows:
        require(isinstance(row, Mapping), "source file row is not an object")
        key = row.get("path")
        require(isinstance(key, str) and key not in hashes, "unexpected or duplicate frozen input")
        path = frozen_path(key, pending_path)
        require(path not in separately_authenticated_briefs, "global/pending brief must be authenticated separately, not listed in source manifest files")
        require(path not in resolved_rows, "duplicate frozen path: " + str(path))
        expected = row.get("sha256")
        require(isinstance(expected, str) and len(expected) == 64, "invalid source hash: " + key)
        actual = sha_file(path)
        require(actual == expected.lower(), "frozen source hash changed: " + key)
        require(row.get("gzip_decoded", False) is False, "gzip decode marker must remain false")
        hashes[canonical_key(path)] = actual
        resolved_rows.add(path)
    require(resolved_rows == expected_paths, "root manifest does not freeze exactly required inputs")
    return source, hashes, source_sha, global_sha, pending_sha, pending_path


def authenticated_bytes(
    hashes: Mapping[str, str],
    pending_path: Path,
    global_sha: str,
    pending_sha: str,
) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for key, expected in hashes.items():
        data = frozen_path(key, pending_path).read_bytes()
        require(sha_bytes(data) == expected, "frozen input changed after entry: " + key)
        raw[key] = data
    global_raw = GLOBAL_BRIEF.read_bytes()
    require(sha_bytes(global_raw) == global_sha, "global continuation brief changed after entry")
    pending_raw = pending_path.read_bytes()
    require(sha_bytes(pending_raw) == pending_sha, "global pending brief changed after entry")
    return raw


def static_new_import_guard(raw: Mapping[str, bytes], relative: str) -> dict[str, Any]:
    """Reject old producer imports without importing or compiling the module."""

    try:
        tree = ast.parse(raw[relative].decode("utf-8"), filename=relative)
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise RunnerInputError("new module source is not parseable: " + relative) from exc
    imports: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    forbidden = [name for name in imports if any(name.startswith(prefix) for prefix in FORBIDDEN_IMPORT_PREFIXES)]
    require(not forbidden, "new module imports forbidden historical code: " + ", ".join(forbidden))
    return {"path": relative, "imports": sorted(set(imports)), "old_producer_imports": False, "static_only": True}


def load_new_module(relative: str, name: str) -> Any:
    path = repository_path(relative)
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "new module loader missing: " + relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def error_record(exc: BaseException) -> dict[str, str]:
    return {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}


def fraction_like(value: Any) -> Fraction | None:
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError):
            return None
    if isinstance(value, Mapping) and value.get("kind") == "exact_fraction":
        try:
            return Fraction(int(value["numerator"]), int(value["denominator"]))
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            return None
    return None


def validate_whole_z_result(value: Any) -> dict[str, Any]:
    require(isinstance(value, Mapping), "new leaf result must be a mapping")
    require(value.get("schema") == SCHEMA, "new leaf result schema mismatch")
    require(value.get("whole_Z_coverage") is True, "new leaf did not certify four-cell whole-Z coverage")
    require(value.get("live_owner_claimed") is False, "new leaf overclaimed live owner")
    require(value.get("ns_closed") is False, "new leaf overclaimed NS closure")
    require(value.get("physical_m15") is False, "new leaf overclaimed physical m15")
    require(value.get("delta_upper_used_as_value") is False, "delta upper bound was treated as a value")
    require(value.get("source_identity_proved") is False, "new leaf overclaimed source identity")
    require(value.get("old_integral_recomputed") is False, "old integral recomputation was claimed")
    require(value.get("old_four_cell_result_reused") is False, "old four-cell result reuse was claimed")
    require(value.get("j") is not None and value.get("delta_interval") is not None, "leaf omitted exact j/delta packet")
    cells = value.get("cells")
    require(isinstance(cells, (list, tuple)) and len(cells) == 4, "leaf must return exactly four cells")
    expected = CELLS
    normalized_domains: list[list[str]] = []
    required_fields = ("G", "G_Z", "integral_error", "anchor_error")
    for index, (cell, (left, right)) in enumerate(zip(cells, expected)):
        require(isinstance(cell, Mapping), f"leaf cell {index} is not an object")
        domain = cell.get("domain")
        require(isinstance(domain, (list, tuple)) and len(domain) == 2, f"leaf cell {index} domain missing")
        got_left, got_right = fraction_like(domain[0]), fraction_like(domain[1])
        require(got_left == left and got_right == right, f"leaf cell {index} domain mismatch")
        for field in required_fields:
            require(field in cell, f"leaf cell {index} missing {field}")
        normalized_domains.append([str(left), str(right)])
    anchor = value.get("anchor")
    require(isinstance(anchor, Mapping), "leaf anchor packet missing")
    require(anchor.get("unique_root_on_domain") is True, "leaf did not retain unique H root premise")
    anchor_domain = anchor.get("domain")
    require(isinstance(anchor_domain, (list, tuple)) and len(anchor_domain) == 2, "leaf anchor domain missing")
    require(fraction_like(anchor_domain[0]) == ANCHOR_DOMAIN[0] and fraction_like(anchor_domain[1]) == ANCHOR_DOMAIN[1], "leaf anchor domain mismatch")
    return {"cell_domains": normalized_domains, "required_fields": list(required_fields), "anchor_domain": [str(ANCHOR_DOMAIN[0]), str(ANCHOR_DOMAIN[1])]}


def consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    pending_path: Path,
    global_sha: str,
    pending_sha: str,
    progress: dict[str, Any],
) -> dict[str, Any]:
    raw = authenticated_bytes(hashes, pending_path, global_sha, pending_sha)
    static_guards = [static_new_import_guard(raw, LEAF), static_new_import_guard(raw, FIXTURE)]
    progress["static_import_guards"] = static_guards

    # The real leaf exposes an owner class and per-cell evaluate/evaluate_cell
    # methods.  Import the new module once, then pass that module to the
    # fixture's sole frozen entry point.
    leaf = load_new_module(LEAF, COMPONENT + "_leaf")
    fixture = load_new_module(FIXTURE, COMPONENT + "_fixture")
    fixture_runner = getattr(fixture, "run_fixture", None)
    require(callable(fixture_runner), "new fixture must export run_fixture(leaf)")
    signature = inspect.signature(fixture_runner)
    parameters = list(signature.parameters.values())
    require(
        len(parameters) == 1
        and parameters[0].kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD),
        "fixture run_fixture signature is outside the frozen run_fixture(leaf) contract",
    )
    try:
        fixture_result = fixture_runner(leaf)
        fixture_api_used = "run_fixture(leaf)"
    except BaseException as exc:
        progress["fixture_calls"] = 1
        progress["fixture_error"] = error_record(exc)
        raise
    progress["fixture_calls"] = 1
    progress["fixture_api_used"] = fixture_api_used
    progress["fixture"] = json_safe(fixture_result)
    require(isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True, "new fixture failed")

    owner_type = getattr(leaf, "ActualCoreGaugeWholeZ", None)
    source_type = getattr(leaf, "GaugeSource", None)
    interval_type = getattr(leaf, "RationalInterval", None)
    require(callable(owner_type) and callable(source_type) and callable(interval_type), "new leaf owner/source API missing")
    source_object = source_type(
        j=J,
        delta=interval_type(DELTA_INTERVAL[0], DELTA_INTERVAL[1]),
        sigma=J / 500,
    )
    owner = owner_type(source=source_object)
    progress["leaf_import_calls"] = 1
    progress["leaf_owner_calls"] = 1
    evaluate_cell = getattr(owner, "evaluate_cell", None)
    evaluate = getattr(owner, "evaluate", None)
    require(callable(evaluate_cell) or callable(evaluate), "new leaf owner has no evaluate/evaluate_cell API")
    raw_cells: list[Mapping[str, Any]] = []
    cell_method = "evaluate_cell" if callable(evaluate_cell) else "evaluate"
    for left, right in CELLS:
        raw_cell = getattr(owner, cell_method)(left, right)
        require(isinstance(raw_cell, Mapping), "new leaf cell result must be a mapping")
        raw_cells.append(raw_cell)
        progress["leaf_cell_evaluation_calls"] = len(raw_cells)
        progress["partial_leaf_cells"] = json_safe(raw_cells)

    anchor_object = getattr(owner, "anchor", None)
    require(anchor_object is not None and callable(getattr(anchor_object, "record", None)), "new leaf anchor certificate missing")
    source_record = source_object.record()
    anchor_record = anchor_object.record()
    normalized_cells: list[dict[str, Any]] = []
    for (left, right), raw_cell in zip(CELLS, raw_cells):
        require(raw_cell.get("schema") == "original_core_gauge_wholez_v1", "new leaf cell schema mismatch")
        require(isinstance(raw_cell.get("G_interval"), Mapping), "new leaf G interval missing")
        require(isinstance(raw_cell.get("G_Z_interval"), Mapping), "new leaf G_Z interval missing")
        budget = raw_cell.get("error_budget")
        require(isinstance(budget, Mapping), "new leaf error budget missing")
        require("integration_truncation" in budget and "anchor_uncertainty" in budget, "new leaf integral/anchor errors missing")
        require(raw_cell.get("source_owner_proved") is False, "new leaf overclaimed source owner")
        require(raw_cell.get("physical_m15") is False and raw_cell.get("closure") is False, "new leaf overclaimed physical closure")
        require(raw_cell.get("delta_center_selected") is False and raw_cell.get("anchor_center_selected") is False, "new leaf selected an enclosure center")
        require(raw_cell.get("midpoint_sample_used") is False, "new leaf used a midpoint sample")
        normalized_cells.append({
            "domain": [left, right],
            "G": raw_cell["G_interval"],
            "G_Z": raw_cell["G_Z_interval"],
            "integral_error": budget["integration_truncation"],
            "anchor_error": budget["anchor_uncertainty"],
            "raw": raw_cell,
        })
    leaf_result = {
        "schema": SCHEMA,
        "j": source_record["j"],
        "delta_interval": source_record["delta_interval"],
        "anchor": {
            **anchor_record,
            "domain": [ANCHOR_DOMAIN[0], ANCHOR_DOMAIN[1]],
            "unique_root_on_domain": anchor_record.get("unique_root_proved_for_every_delta_in_interval") is True,
        },
        "cells": normalized_cells,
        "whole_Z_coverage": True,
        "live_owner_claimed": False,
        "source_identity_proved": False,
        "delta_upper_used_as_value": False,
        "physical_m15": False,
        "ns_closed": False,
        "old_integral_recomputed": False,
        "old_four_cell_result_reused": False,
        "owner_class": "ActualCoreGaugeWholeZ",
        "cell_method": cell_method,
        "source_record": source_record,
        "anchor_record": anchor_record,
    }
    validation = validate_whole_z_result(leaf_result)
    return {
        "schema": SCHEMA,
        "source_manifest_schema": MANIFEST_SCHEMA,
        "source_manifest": source,
        "input_hashes": dict(hashes),
        "fixture_calls": 1,
        "leaf_calls": 1,
        "leaf_import_calls": 1,
        "leaf_owner_calls": 1,
        "leaf_cell_evaluation_calls": len(raw_cells),
        "fixture_api_used": fixture_api_used,
        "static_import_guards": static_guards,
        "fixture": json_safe(fixture_result),
        "leaf_validation": validation,
        "leaf_result": json_safe(leaf_result),
        "parameters": {
            "j": [J.numerator, J.denominator],
            "delta_interval": [
                [DELTA_INTERVAL[0].numerator, DELTA_INTERVAL[0].denominator],
                [DELTA_INTERVAL[1].numerator, DELTA_INTERVAL[1].denominator],
            ],
            "anchor_domain": [[ANCHOR_DOMAIN[0].numerator, ANCHOR_DOMAIN[0].denominator], [0, 1]],
            "cells": [[[left.numerator, left.denominator], [right.numerator, right.denominator]] for left, right in CELLS],
            "arithmetic": "exact Fraction with leaf-directed dyadic outward rounding",
        },
        "coverage": {
            "whole_Z": True,
            "cells": 4,
            "domain": "[-1,1]",
            "ordinary_Z_orders": [0, 1],
            "function_points_materialized": False,
        },
        "scope": {
            "live_owner_claimed": False,
            "source_identity_proved": False,
            "delta_upper_used_as_value": False,
            "physical_m15": False,
            "ns_closed": False,
            "old_producer_called": False,
            "old_fixture_called": False,
            "old_integral_recomputed": False,
            "old_four_cell_result_reused": False,
        },
    }


def main() -> int:
    output: dict[str, Any] = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "phase": PHASE,
        "status": "failed",
        "fixture_calls": 0,
        "leaf_calls": 0,
        "briefs_separately_authenticated": False,
        "whole_Z_coverage": False,
        "live_owner_claimed": False,
        "source_identity_proved": False,
        "delta_upper_used_as_value": False,
        "physical_m15": False,
        "ns_closed": False,
    }
    entry_written = False
    try:
        # The root agent must freeze source.json and the global pending brief;
        # this runner never creates a substitute local brief.
        source, hashes, source_sha, global_sha, pending_sha, pending_path = preflight_frozen()
        entry_sha = save_exclusive(
            ENTRY,
            {
                "schema": COMPONENT + "_entry_v1",
                "component": COMPONENT,
                "phase": PHASE,
                "source_sha256": source_sha,
                "global_brief_sha256": global_sha,
                "pending_brief_sha256": pending_sha,
                "input_hashes": hashes,
                "entered": True,
                "fixture_calls": 0,
                "leaf_calls": 0,
                "no_gzip_decode": True,
                "briefs_separately_authenticated": True,
                "old_producer_called": False,
            },
        )
        entry_written = True
        output.update({"source_sha256": source_sha, "global_brief_sha256": global_sha, "pending_brief_sha256": pending_sha, "entry_sha256": entry_sha, "input_hashes": hashes, "briefs_separately_authenticated": True})
        output.update(consume(source, hashes, pending_path, global_sha, pending_sha, output))
        output["status"] = "returned"
    except BaseException as exc:
        output["error"] = error_record(exc)
    if entry_written:
        result_sha = save_exclusive(RESULT, output)
        print(json.dumps({"component": COMPONENT, "status": output["status"], "entry_sha256": output.get("entry_sha256"), "result_sha256": result_sha}, sort_keys=True))
        return 0 if output["status"] == "returned" else 1
    # A startup/freeze failure happened before the one-entry boundary.  Keep
    # its traceback separately, without fabricating an entered component.
    startup = DIRECTORY / "startup_failure.json"
    if not startup.exists():
        save_exclusive(startup, output)
    print(json.dumps({"component": COMPONENT, "status": "startup_failed", "startup_failure": startup.as_posix()}, sort_keys=True))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
