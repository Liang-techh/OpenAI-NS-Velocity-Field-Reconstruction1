"""Enter the F-gap radius-spine source capsule exactly once.

Every raw file named by the pending source manifest is authenticated before
the exclusive entry receipt; the two gzip archives are hashed as bytes only
at that point.  After entry, this runner imports only the new radius-spine
leaf and invokes its two public interfaces once: the synthetic fixture and
the bounded source-spine builder.

The source builder is allowed to read three distinct indexed positions in
each of the two archives (six archive-row lookups total) and decode the two
authenticated gzip archives.  That is a bounded raw-row read, not a claim
that no graph rows were accessed.  The runner records both counts in the
success receipt and rejects an unbounded traversal, repeated graph
computation, source-value evaluation, or integration.

Both success and post-entry failure paths persist a complete JSON receipt
before printing a short summary.  Failure receipts retain any JSON-safe leaf
result already returned by either interface.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "f_gap_radius_spine_binding_v1"
PHASE = f"{COMPONENT}_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "f_gap_radius_spine_binding.py"

ARCHIVE_RP_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
ARCHIVE_N1_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
)
HISTORICAL_DEFERRED_REL = (
    "tmp/historical_frame_gap_deferred_roles_v1/result.json"
)
HISTORICAL_REWRITE_REL = (
    "tmp/historical_frame_gap_saved_rewrite_lookup_v1/result.json"
)
OWN_FIELDS_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py"
)
RP_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py"
)
TRANSPORT_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_Rc_function_transport.py"
)

REQUIRED_INPUTS = {
    "src/openai_ns_reconstruction/f_gap_radius_spine_binding.py",
    "scripts/run_f_gap_radius_spine_binding_once.py",
    HISTORICAL_DEFERRED_REL,
    HISTORICAL_REWRITE_REL,
    ARCHIVE_RP_REL,
    ARCHIVE_N1_REL,
    OWN_FIELDS_REL,
    RP_SOURCE_REL,
    TRANSPORT_REL,
}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


def sha256(path: Path) -> str:
    """Hash raw bytes in bounded chunks; do not decode archives here."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _json_bytes(value: Any) -> bytes:
    """Canonical JSON bytes, rejecting NaN and non-JSON leaf values."""

    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def write_once(path: Path, value: Mapping[str, Any]) -> str:
    """Persist one receipt exclusively and return its persisted-byte SHA."""

    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _repository_path(raw_path: Any, label: str = "manifest path") -> tuple[str, Path]:
    require(isinstance(raw_path, str) and bool(raw_path), label + " is missing")
    relative = Path(raw_path)
    require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " is missing or escapes repository: " + raw_path,
    )
    return key, path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate every direct raw path/hash record before entry.

    Checkpoint writers use either ``files: [{path, sha256}, ...]`` or a
    direct hash mapping.  Both layouts are accepted, but no nested source
    graph is walked while collecting the manifest.
    """

    records: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            require(
                isinstance(item, Mapping),
                f"source.files[{index}] is not an object",
            )
            records.append(
                (item.get("path"), item.get("sha256"), f"source.files[{index}]")
            )
    for key in ("input_hashes", "sha256", "input_sha256", "raw_input_sha256", "code_sha256"):
        mapping = source.get(key)
        if isinstance(mapping, Mapping):
            for raw_path, expected in mapping.items():
                records.append((raw_path, expected, f"source.{key}"))
    require(bool(records), "source manifest has no direct raw file hashes")

    hashes: dict[str, str] = {}
    for raw_path, expected, label in records:
        key, path = _repository_path(raw_path, label + " path")
        require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        digest = expected.lower()
        require(key not in hashes, "duplicate manifest path: " + key)
        actual = sha256(path)
        require(actual.lower() == digest, "raw SHA mismatch: " + key)
        hashes[key] = actual

    required = {item.lower() for item in REQUIRED_INPUTS}
    present = {item.lower() for item in hashes}
    missing = sorted(item for item in required if item not in present)
    require(
        not missing,
        "manifest omits required frozen input(s): " + ", ".join(missing),
    )
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    """Authenticate the pending brief's pointer to this source manifest."""

    checkpoint: Any = brief.get(f"{COMPONENT}_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    if isinstance(checkpoint, Mapping):
        raw_path = checkpoint.get("path")
        expected = checkpoint.get("sha256")
    else:
        raw_path = brief.get("source_checkpoint_path")
        expected = brief.get("source_checkpoint_sha256")
    require(isinstance(raw_path, str), "pending brief source checkpoint path missing")
    checkpoint_path = Path(raw_path)
    if not checkpoint_path.is_absolute():
        checkpoint_path = ROOT / checkpoint_path
    require(
        checkpoint_path.resolve() == SOURCE.resolve(),
        "pending brief source checkpoint points to another file",
    )
    actual = sha256(SOURCE)
    require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and actual.lower() == expected.lower(),
        "pending brief source checkpoint SHA-256 mismatch",
    )
    return actual


def _count_for_keys(marker: Any, keys: tuple[str, ...]) -> int | None:
    """Read one explicit count from a scalar/list/mapping marker."""

    if isinstance(marker, bool):
        return None
    if isinstance(marker, int):
        return marker
    if isinstance(marker, (list, tuple)):
        return len(marker)
    if isinstance(marker, Mapping):
        for key in keys:
            candidate = marker.get(key)
            if isinstance(candidate, int) and not isinstance(candidate, bool):
                return candidate
            if isinstance(candidate, (list, tuple, Mapping)):
                return len(candidate)
    return None


def _bounded_counts(value: Mapping[str, Any]) -> tuple[int, int, Any] | None:
    """Return distinct-index and per-archive lookup counts from scope data."""

    names = (
        "bounded_source_row_read",
        "bounded_source_rows_read",
        "bounded_source_index_read",
        "bounded_archive_row_read",
        "source_row_read",
        "source_rows_accessed",
        "archive_index_rows",
        "archive_index_nodes",
        "graph_rows_accessed",
        "source_graph_rows_accessed",
    )
    containers = (value,)
    for container_name in (
        "scope",
        "scope_flags",
        "source_scope",
        "bounded_scope",
        "archive_scope",
    ):
        nested = value.get(container_name)
        if isinstance(nested, Mapping):
            containers += (nested,)

    distinct_keys = (
        "distinct_archive_indices",
        "distinct_index_count",
        "distinct_indices",
        "archive_index_count",
        "index_count",
        "indices",
        "indexes",
    )
    lookup_keys = (
        "archive_row_lookups",
        "archive_row_lookup_count",
        "row_lookups",
        "lookup_count",
        "lookups",
        "total_row_count",
        "row_count",
    )
    for container in containers:
        for name in names:
            if name not in container:
                continue
            marker = container[name]
            if not isinstance(marker, Mapping):
                continue
            distinct = _count_for_keys(marker, distinct_keys)
            lookups = _count_for_keys(marker, lookup_keys)
            if distinct == 3 and lookups == 6:
                return distinct, lookups, marker
        # Some leaves expose the two counters as sibling fields rather than
        # nesting them under bounded_source_row_read.
        distinct = _count_for_keys(container, distinct_keys)
        lookups = _count_for_keys(container, lookup_keys)
        if distinct == 3 and lookups == 6:
            return distinct, lookups, container
    return None


def _bounded_manifest_guard(source: Mapping[str, Any]) -> None:
    """Require three distinct indices and six archive-row lookups."""

    counts = _bounded_counts(source)
    require(
        counts is not None,
        "source manifest must record distinct_indices=3 and archive_row_lookups=6",
    )


def preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    """Read pending metadata and authenticate all declared raw bytes."""

    require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    require(BRIEF.is_file(), "pending brief is missing")
    require(SOURCE.is_file(), "source manifest is missing")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    require(isinstance(brief, Mapping), "pending brief is not an object")
    require(isinstance(source, Mapping), "source manifest is not an object")
    require(
        brief.get("component") == COMPONENT and brief.get("phase") == PHASE,
        "pending brief component or phase mismatch",
    )
    require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal flags mismatch",
    )
    source_sha = _checkpoint_sha(brief)
    require(source.get("component") == COMPONENT, "source manifest component mismatch")
    require(source.get("one_entry_only") is True, "source manifest is not one-entry-only")
    require(
        source.get("all_raw_inputs_authenticated_before_decompression") is True,
        "raw-before-decompression guard is missing",
    )
    require(
        source.get("no_source_graph_traversal") is True,
        "source graph traversal guard is missing",
    )
    _bounded_manifest_guard(source)
    for name in (
        "no_duplicate_graph_computation",
        "no_repeated_graph_computation",
        "no_old_leaf_import",
        "no_old_producer_import",
    ):
        if name in source:
            require(source[name] is True, "source manifest guard is false: " + name)
    hashes = _manifest_hashes(source)
    require(LEAF.is_file(), "new radius-spine leaf is missing")
    return dict(source), hashes, source_sha


def _load_new_module() -> Any:
    """Import only the new leaf, after the exclusive entry receipt."""

    spec = importlib.util.spec_from_file_location(
        "f_gap_radius_spine_binding_new_leaf", LEAF
    )
    require(spec is not None and spec.loader is not None, "new radius-spine leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _require_false_group(
    value: Mapping[str, Any],
    names: tuple[str, ...],
    label: str,
    *,
    required: bool = True,
) -> None:
    present = [name for name in names if name in value]
    if required:
        require(bool(present), label + " scope flag is missing")
    for name in present:
        require(value[name] is False, label + " scope flag is true: " + name)


def _require_scope_flags(value: Any, label: str) -> None:
    require(isinstance(value, Mapping), label + " did not return an object")
    _require_false_group(
        value,
        ("source_values_computed", "source_value_computation", "actual_source_values_computed"),
        label + " source-value",
    )
    _require_false_group(
        value,
        ("integrals_evaluated", "integrals_computed", "quadrature_evaluated", "integral_calls"),
        label + " integral",
    )
    _require_false_group(
        value,
        ("numeric_use_admissible", "admissible_for_numeric_use"),
        label + " numeric-use",
    )
    _require_false_group(value, ("certified",), label + " certification")
    _require_false_group(
        value,
        ("source_graph_traversed", "graph_traversed"),
        label + " graph-traversal",
    )
    _require_false_group(
        value,
        ("duplicate_graph_computation", "repeated_graph_computation", "graph_recomputed"),
        label + " duplicate-computation",
        required=False,
    )
    _require_false_group(
        value,
        ("old_leaf_imported", "old_producer_called", "ancestor_math_called", "ancestor_replayed"),
        label + " old-code",
        required=False,
    )


def _require_fixture(value: Any) -> None:
    require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    markers = [name for name in ("passed", "checks_passed", "fixture_passed") if name in value]
    require(bool(markers), "synthetic fixture omitted a pass marker")
    require(all(value[name] is True for name in markers), "synthetic fixture did not pass")
    _require_scope_flags(value, "synthetic fixture")
    _json_bytes(value)


def _count_from_marker(marker: Any) -> int | None:
    """Extract a count from a JSON-safe row/archive evidence marker."""

    if isinstance(marker, bool):
        return None
    if isinstance(marker, int):
        return marker
    if isinstance(marker, (list, tuple)):
        return len(marker)
    if isinstance(marker, Mapping):
        for key in (
            "count",
            "row_count",
            "index_count",
            "access_count",
            "archive_count",
            "decompressed_count",
            "parsed_count",
        ):
            candidate = marker.get(key)
            if isinstance(candidate, int) and not isinstance(candidate, bool):
                return candidate
            if isinstance(candidate, (list, tuple)):
                return len(candidate)
        for key in ("rows", "indices", "indexes", "index_nodes", "archives", "paths"):
            candidate = marker.get(key)
            if isinstance(candidate, (list, tuple, Mapping)):
                return len(candidate)
    return None


def _find_evidence(
    value: Mapping[str, Any],
    names: tuple[str, ...],
    *,
    expected: int,
    label: str,
) -> dict[str, Any]:
    """Find direct or one-level nested evidence for an exact bounded count."""

    candidates: list[tuple[str, Any]] = []
    for name in names:
        if name in value:
            candidates.append((name, value[name]))
    for container_name in (
        "scope",
        "scope_flags",
        "source_scope",
        "bounded_scope",
        "archive_scope",
    ):
        nested = value.get(container_name)
        if isinstance(nested, Mapping):
            for name in names:
                if name in nested:
                    candidates.append((f"{container_name}.{name}", nested[name]))
    for name, marker in candidates:
        count = _count_from_marker(marker)
        if count == expected:
            return {"field": name, "count": count, "evidence": marker}
    require(False, f"{label} evidence must record exactly {expected}")
    raise AssertionError("unreachable")


def _require_source_spine(value: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    require(isinstance(value, Mapping), "bounded source spine did not return an object")
    _require_scope_flags(value, "bounded source spine")
    counts = _bounded_counts(value)
    require(
        counts is not None,
        "bounded source-row evidence must record distinct_indices=3 "
        "and archive_row_lookups=6",
    )
    distinct_indices, archive_row_lookups, marker = counts
    row_evidence = {
        "distinct_archive_indices": distinct_indices,
        "archive_row_lookups": archive_row_lookups,
        "evidence": marker,
    }
    archive_evidence = _find_evidence(
        value,
        (
            "gzip_archives_parsed",
            "gzip_archives_decompressed",
            "archives_parsed",
            "archive_parse_count",
            "archive_decompression_count",
            "decompressed_archives",
            "parsed_archives",
            "archive_paths",
        ),
        expected=2,
        label="gzip archive parse",
    )
    # A bounded row read is intentionally not represented as zero graph-row
    # access.  If a leaf emits such a contradictory marker, reject it.
    for name in (
        "no_graph_rows_accessed",
        "graph_rows_accessed_zero",
        "source_graph_rows_accessed_zero",
    ):
        if name in value:
            require(value[name] is False, "source spine contradicts bounded row evidence: " + name)
    _json_bytes(value)
    return row_evidence, archive_evidence


def _recordable(value: Any) -> Any:
    """Retain a returned leaf packet in failure output only if JSON-safe."""

    if value is None:
        return None
    try:
        _json_bytes(value)
    except Exception:
        return None
    return value


def main() -> int:
    source, hashes, source_sha = preflight()
    entry_sha = write_once(
        ENTRY,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "consumer_calls_started": 0,
            "verification_invocations": 0,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "graph_row_scope": None,
            "distinct_archive_indices": None,
            "archive_row_lookups": 0,
            "bounded_source_row_read": None,
            "gzip_archives_parsed": 0,
            "integrals_evaluated": False,
            "duplicate_graph_computation": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
        },
    )

    calls: list[str] = []
    invocation_count = 0
    fixture_result: Any = None
    spine_result: Any = None
    row_evidence: Any = None
    archive_evidence: Any = None
    try:
        module = _load_new_module()
        fixture = getattr(module, "verify_synthetic_fixture", None)
        require(callable(fixture), "new radius-spine leaf has no verify_synthetic_fixture")
        calls.append("verify_synthetic_fixture")
        invocation_count += 1
        fixture_result = fixture()
        _require_fixture(fixture_result)

        builder = getattr(module, "build_radius_spine_binding", None)
        require(callable(builder), "new radius-spine leaf has no build_radius_spine_binding")
        calls.append("build_radius_spine_binding")
        invocation_count += 1
        spine_result = builder(ROOT, source, hashes)
        row_evidence, archive_evidence = _require_source_spine(spine_result)

        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": calls,
            "verification_invocations": invocation_count,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "graph_row_scope": {
                "distinct_archive_indices": row_evidence["distinct_archive_indices"],
                "archive_row_lookups": row_evidence["archive_row_lookups"],
            },
            "distinct_archive_indices": row_evidence["distinct_archive_indices"],
            "archive_row_lookups": row_evidence["archive_row_lookups"],
            "bounded_source_row_read": row_evidence,
            "gzip_archives_parsed": 2,
            "archive_parse_evidence": archive_evidence,
            "integrals_evaluated": False,
            "duplicate_graph_computation": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
            "result": {
                "synthetic_fixture": fixture_result,
                "radius_spine_binding": spine_result,
            },
            "error": None,
        }
        # Validate the complete success packet before its exclusive write.
        _json_bytes(output)
    except Exception as exc:  # complete post-entry failure artifact is required
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": calls,
            "verification_invocations": invocation_count,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "graph_row_scope": (
                {
                    "distinct_archive_indices": row_evidence["distinct_archive_indices"],
                    "archive_row_lookups": row_evidence["archive_row_lookups"],
                }
                if row_evidence is not None
                else None
            ),
            "distinct_archive_indices": (
                row_evidence["distinct_archive_indices"] if row_evidence is not None else None
            ),
            "archive_row_lookups": (
                row_evidence["archive_row_lookups"] if row_evidence is not None else 0
            ),
            "bounded_source_row_read": _recordable(row_evidence),
            "gzip_archives_parsed": 2 if archive_evidence is not None else 0,
            "archive_parse_evidence": _recordable(archive_evidence),
            "integrals_evaluated": False,
            "duplicate_graph_computation": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
            "fixture_result": _recordable(fixture_result),
            "radius_spine_result": _recordable(spine_result),
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        }

    result_sha = write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": calls,
                "verification_invocations": invocation_count,
                "distinct_archive_indices": output["distinct_archive_indices"],
                "archive_row_lookups": output["archive_row_lookups"],
                "gzip_archives_parsed": output["gzip_archives_parsed"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
