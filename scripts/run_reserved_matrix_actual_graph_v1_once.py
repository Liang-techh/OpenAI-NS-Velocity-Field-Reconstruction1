"""Run the reserved actual graph and same-graph C1 action once.

The root agent must freeze ``tmp/reserved_matrix_actual_graph_v1/source.json``
before invoking this runner.  The manifest pins the new graph consumer, the
injected C1 action, this runner, and exactly two gzip archives marked
``reserved`` and ``full``.  Every raw byte is authenticated before the
exclusive entry receipt is written.  Only after entry are the new modules
imported and the archives decoded.

This runner deliberately does not import an old producer, call a fixture, or
evaluate a graph handle.  The result records opaque handles and bounded
structural comparisons with the decoded archive rows/coefficient containers;
it does not claim physical closure, inverse verification, or a numerical
value.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib
import json
import re
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
NAME = "reserved_matrix_actual_graph_v1"
DIRECTORY = ROOT / "tmp" / NAME
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

ACTUAL_GRAPH = "src/openai_ns_reconstruction/reserved_matrix_actual_graph_v1.py"
C1_ACTION = "src/openai_ns_reconstruction/reserved_matrix_c1_action_v1.py"
RUNNER = "scripts/run_reserved_matrix_actual_graph_v1_once.py"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

_ARCHIVE_ROLES = ("reserved", "full")
_MISSING = object()
_SCALAR_TYPES = (str, int, float, bool, type(None))
_PAIR_FIELDS = ("value", "Z")
_HANDLE_ID_FIELDS = ("node", "function_id", "handle_id", "source_node")
_ROW_KEYS = ("exact_graph_nodes", "graph_nodes", "nodes", "rows")
_COEFFICIENT_KEY_ALIASES = {
    "m",
    "mz",
    "minverse",
    "matrix",
    "matrixz",
    "matrixinverse",
    "y",
    "yz",
    "coefficients",
    "coefficientmatrix",
    "coefficientrows",
    "matrixcoefficients",
}


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is malformed or changed."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON value is not admitted: " + token)


def _load_json(raw: bytes, label: str) -> Any:
    try:
        return json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite UTF-8 JSON") from exc


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


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    try:
        with path.open("xb") as stream:
            stream.write(raw)
            stream.flush()
    except FileExistsError as exc:
        raise RunnerInputError("one-entry file already exists: " + str(path)) from exc
    return _sha(raw)


def _repository_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(), label + " is missing")
    return key, path


def _expected_sha(row: Mapping[str, Any], label: str) -> str:
    raw = row.get("raw_sha256", row.get("sha256"))
    if row.get("raw_sha256") is not None and row.get("sha256") is not None:
        _require(
            str(row["raw_sha256"]).lower() == str(row["sha256"]).lower(),
            label + " SHA fields disagree",
        )
    _require(isinstance(raw, str) and SHA256_RE.fullmatch(raw) is not None,
             label + " has invalid SHA-256")
    return raw.lower()


def _normalise_role(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    token = re.sub(r"[^a-z0-9]", "", value.lower())
    if token in {"reserved", "reservedarchive", "reservedgraph"}:
        return "reserved"
    if token in {"full", "fullarchive", "fullgraph"}:
        return "full"
    return None


def _archive_declarations(source: Mapping[str, Any]) -> dict[str, tuple[str, str, str]]:
    """Return role -> (repository path, expected raw SHA, declaration label).

    The preferred frozen form is a ``role`` field on each archive row in
    ``source.files``.  A small equivalent ``source.archives`` mapping/list is
    accepted so the root agent can keep archive metadata separate from the
    ordinary frozen source files.  Roles remain explicit; a filename guess is
    never used for a real archive.
    """

    found: dict[str, tuple[str, str, str]] = {}

    def add(role_value: Any, row: Mapping[str, Any], label: str) -> None:
        role = _normalise_role(role_value)
        if role is None:
            return
        key, _ = _repository_path(row.get("path"), label)
        expected = _expected_sha(row, label)
        _require(key.lower().endswith(".gz"), label + " must be a gzip archive")
        _require(
            row.get("raw_sha256_before_decode") is True
            or row.get("raw_bytes_authenticated_before_decode") is True,
            label + " is not marked raw-authenticated before decode",
        )
        previous = found.get(role)
        item = (key, expected, label)
        if previous is not None:
            _require(previous[:2] == item[:2], "conflicting " + role + " archive declarations")
        else:
            found[role] = item

    files = source.get("files")
    _require(isinstance(files, list), "source.files is missing")
    for index, row in enumerate(files):
        _require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        role = row.get("role", row.get("archive_role", row.get("archive")))
        add(role, row, f"source.files[{index}]")

    archives = source.get("archives")
    if isinstance(archives, Mapping):
        for raw_role, raw_row in archives.items():
            _require(isinstance(raw_row, Mapping), f"source.archives[{raw_role!r}] is not an object")
            add(raw_role, raw_row, f"source.archives[{raw_role!r}]")
    elif archives is not None:
        _require(isinstance(archives, list), "source.archives must be a mapping or list")
        for index, row in enumerate(archives):
            _require(isinstance(row, Mapping), f"source.archives[{index}] is not an object")
            add(row.get("role", row.get("archive_role")), row, f"source.archives[{index}]")

    for raw_role in _ARCHIVE_ROLES:
        top_level = source.get(raw_role + "_archive")
        if isinstance(top_level, Mapping):
            add(raw_role, top_level, "source." + raw_role + "_archive")
        elif top_level is not None:
            raise RunnerInputError("source." + raw_role + "_archive must be an object")

    _require(set(found) == set(_ARCHIVE_ROLES),
             "source must explicitly declare exactly reserved and full archives")
    _require(found["reserved"][0] != found["full"][0],
             "reserved and full archives must be distinct files")
    return found


def _manifest_hashes(
    source: Mapping[str, Any],
    archives: Mapping[str, tuple[str, str, str]],
) -> dict[str, str]:
    """Authenticate every frozen file and both archive raw byte streams."""

    rows = source.get("files")
    _require(isinstance(rows, list), "source.files is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        label = f"source.files[{index}]"
        _require(isinstance(row, Mapping), label + " is not an object")
        key, path = _repository_path(row.get("path"), label)
        _require(key not in hashes, "duplicate source file: " + key)
        expected = _expected_sha(row, label)
        actual = _sha_file(path)
        _require(actual.lower() == expected, "raw SHA-256 mismatch: " + key)
        if key.lower().endswith(".gz"):
            _require(
                row.get("raw_sha256_before_decode") is True
                or row.get("raw_bytes_authenticated_before_decode") is True,
                "gzip input is not marked raw-authenticated: " + key,
            )
        hashes[key] = actual

    for role, (key, expected, label) in archives.items():
        path = (ROOT / key).resolve()
        _require(path.is_file(), label + " archive is missing")
        actual = _sha_file(path)
        _require(actual.lower() == expected, "raw SHA-256 mismatch: " + role + " archive")
        prior = hashes.get(key)
        if prior is not None:
            _require(prior.lower() == actual.lower(), "archive hash disagrees with source.files: " + key)
        hashes[key] = actual

    required = {ACTUAL_GRAPH, C1_ACTION, RUNNER}
    _require(required.issubset(hashes), "source manifest is missing a frozen new input")
    for role, (key, _, _) in archives.items():
        _require(key in hashes, "archive was not authenticated: " + role)
    return hashes


def _preflight() -> tuple[dict[str, Any], dict[str, str], dict[str, tuple[str, str, str]], str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(SOURCE.is_file(), "source manifest is missing")
    source_raw = SOURCE.read_bytes()
    source = _load_json(source_raw, "source manifest")
    _require(isinstance(source, Mapping), "source manifest must be an object")
    _require(source.get("component") in (None, NAME), "source manifest component mismatch")
    _require(source.get("one_entry_only") in (None, True), "source one-entry guard is missing")
    archives = _archive_declarations(source)
    hashes = _manifest_hashes(source, archives)
    return dict(source), hashes, archives, _sha(source_raw)


def _recheck_raw(hashes: Mapping[str, str], key: str) -> bytes:
    _require(key in hashes, "unhashed input: " + key)
    raw = (ROOT / key).read_bytes()
    _require(_sha(raw).lower() == hashes[key].lower(), "input changed after entry: " + key)
    return raw


def _import_new(path_key: str, hashes: Mapping[str, str]) -> Any:
    path = (ROOT / path_key).resolve()
    _require(_sha_file(path).lower() == hashes[path_key].lower(),
             "new input changed after entry: " + path_key)
    module_name = "openai_ns_reconstruction." + Path(path_key).stem
    _require(module_name not in sys.modules, "new module already imported: " + module_name)
    source_root = str((ROOT / "src").resolve())
    if source_root not in sys.path:
        sys.path.insert(0, source_root)
    module = importlib.import_module(module_name)
    module_path = _safe_get(module, "__file__")
    _require(
        isinstance(module_path, str) and Path(module_path).resolve() == path,
        "imported module path does not match frozen input: " + path_key,
    )
    return module


def _safe_get(value: Any, attribute: str) -> Any:
    try:
        return getattr(value, attribute)
    except (AttributeError, TypeError):
        return _MISSING


def _is_pair(value: Any) -> bool:
    return _safe_get(value, "value") is not _MISSING and _safe_get(value, "Z") is not _MISSING


def _primitive_descriptor(value: Any) -> Any:
    if type(value) in _SCALAR_TYPES:
        return {"type": type(value).__name__}
    return {"type": type(value).__name__}


def _handle_record(handle: Any) -> dict[str, Any]:
    """Serialize an opaque graph handle without reading its source value."""

    record: dict[str, Any] = {
        "kind": "source_handle",
        "type": type(handle).__module__ + "." + type(handle).__qualname__,
    }
    for field in _HANDLE_ID_FIELDS:
        candidate = _safe_get(handle, field)
        if candidate is not _MISSING and type(candidate) in _SCALAR_TYPES:
            record[field] = candidate
    owner = _safe_get(handle, "graph")
    if owner is _MISSING:
        owner = _safe_get(handle, "source_graph")
    if owner is not _MISSING and owner is not None:
        record["owner"] = {
            "type": type(owner).__module__ + "." + type(owner).__qualname__,
            "identity": id(owner),
        }
    return record


def _opaque_record(value: Any, *, depth: int = 0) -> Any:
    """Bounded JSON representation of handles/pairs and their containers."""

    if depth > 8:
        return {"kind": "depth_limit", "type": type(value).__name__}
    if _is_pair(value):
        return {
            "kind": "c1_pair",
            "value": _opaque_record(_safe_get(value, "value"), depth=depth + 1),
            "Z": _opaque_record(_safe_get(value, "Z"), depth=depth + 1),
        }
    if isinstance(value, Mapping):
        return {
            "kind": "mapping",
            "type": type(value).__name__,
            "keys": [str(key) for key in list(value.keys())[:64]],
            "size": len(value),
        }
    if isinstance(value, (str, bytes, bytearray)) or value is None:
        return _primitive_descriptor(value)
    if isinstance(value, Sequence):
        return {
            "kind": "sequence",
            "type": type(value).__name__,
            "size": len(value),
            "items": [_opaque_record(item, depth=depth + 1) for item in list(value)[:16]],
        }
    return _handle_record(value)


def _report_record(value: Any, fields: Sequence[str]) -> dict[str, Any]:
    """Keep scalar provenance flags from a report without reading handles."""

    result: dict[str, Any] = {
        "type": type(value).__module__ + "." + type(value).__qualname__,
    }
    for field in fields:
        candidate = _safe_get(value, field)
        if candidate is _MISSING:
            continue
        if isinstance(candidate, Mapping):
            result[field] = {str(key): candidate[key] for key in list(candidate.keys())[:64]}
        elif isinstance(candidate, (str, int, float, bool)) or candidate is None:
            result[field] = candidate
        elif isinstance(candidate, Sequence) and not isinstance(candidate, (bytes, bytearray, str)):
            result[field] = [item for item in list(candidate)[:64]]
        else:
            result[field] = {"type": type(candidate).__name__}
    return result


def _source_id_literal(value: Any, *, depth: int = 0) -> Any:
    """Copy pinned integer/string source-id tables without touching handles."""

    if depth > 8:
        return {"type": type(value).__name__, "depth_limit": True}
    if type(value) in (str, int, bool, float) or value is None:
        return value
    if isinstance(value, Mapping):
        return {
            str(key): _source_id_literal(child, depth=depth + 1)
            for key, child in value.items()
        }
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray, str)):
        return [_source_id_literal(child, depth=depth + 1) for child in value]
    return {"type": type(value).__module__ + "." + type(value).__qualname__}


def _shape(value: Any, *, depth: int = 0) -> Any:
    """Shape only; scalar contents are intentionally excluded."""

    if depth > 4:
        return {"kind": "depth_limit", "type": type(value).__name__}
    if isinstance(value, Mapping):
        keys = [str(key) for key in value.keys()]
        result: dict[str, Any] = {
            "kind": "mapping",
            "type": type(value).__name__,
            "size": len(value),
            "keys": sorted(keys)[:64],
        }
        selected = {}
        for key in list(value.keys())[:32]:
            selected[str(key)] = _shape(value[key], depth=depth + 1)
        result["selected"] = selected
        return result
    if isinstance(value, (str, bytes, bytearray)) or value is None:
        return {"kind": "scalar", "type": type(value).__name__}
    if isinstance(value, Sequence):
        items = list(value)
        return {
            "kind": "sequence",
            "type": type(value).__name__,
            "size": len(items),
            "first": _shape(items[0], depth=depth + 1) if items else None,
            "last": _shape(items[-1], depth=depth + 1) if items else None,
        }
    return {"kind": "opaque", "type": type(value).__module__ + "." + type(value).__qualname__}


def _shape_equal(left: Any, right: Any) -> bool:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        # Ignore concrete container class names and selected scalar field names;
        # compare only structural kind, sizes, and nested shape.
        if left.get("kind") != right.get("kind"):
            return False
        if left.get("size") != right.get("size"):
            return False
        for key in ("first", "last"):
            if key in left or key in right:
                if not _shape_equal(left.get(key), right.get(key)):
                    return False
        if "selected" in left or "selected" in right:
            lsel = left.get("selected", {})
            rsel = right.get("selected", {})
            if set(lsel) != set(rsel):
                return False
            if any(not _shape_equal(lsel[key], rsel[key]) for key in lsel):
                return False
        return True
    return left == right


def _rows_summary(rows: Any) -> dict[str, Any] | None:
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes, bytearray)):
        return None
    count = len(rows)

    def row_shape(row: Any) -> Any:
        if isinstance(row, Mapping):
            operation = row.get("operation")
            return {
                "kind": "mapping",
                "keys": sorted(str(key) for key in row.keys()),
                "operation_type": type(operation).__name__,
            }
        return {"kind": "other", "type": type(row).__name__}

    first = row_shape(rows[0]) if count else None
    last = row_shape(rows[-1]) if count else None
    prefix_rows = list(rows[:32])
    suffix_rows = list(rows[-32:]) if count > 32 else prefix_rows
    operations: list[str] = []
    suffix_operations: list[str] = []
    for row in prefix_rows:
        if isinstance(row, Mapping) and isinstance(row.get("operation"), str):
            operations.append(row["operation"])
    for row in suffix_rows:
        if isinstance(row, Mapping) and isinstance(row.get("operation"), str):
            suffix_operations.append(row["operation"])
    return {
        "count": count,
        "first_row": first,
        "last_row": last,
        "operation_prefix": operations,
        "operation_suffix": suffix_operations,
        "complete_scan": False,
    }


def _first_archive_rows(archive: Mapping[str, Any]) -> tuple[str | None, Any]:
    for key in _ROW_KEYS:
        value = archive.get(key)
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            return key, value
    graph = archive.get("graph")
    if isinstance(graph, Mapping):
        for key in _ROW_KEYS:
            value = graph.get(key)
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
                return "graph." + key, value
    return None, None


def _graph_rows(graph: Any) -> tuple[str | None, Any]:
    for key in ("nodes", "rows", "graph_nodes"):
        value = _safe_get(graph, key)
        if value is not _MISSING and isinstance(value, Sequence):
            return key, value
    return None, None


def _normalise_key(key: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(key).lower())


def _coefficient_candidates(archive: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Collect shallow coefficient containers without walking graph rows."""

    found: list[dict[str, Any]] = []

    def visit(mapping: Mapping[str, Any], path: str, depth: int) -> None:
        if depth > 3:
            return
        for key, value in mapping.items():
            normalized = _normalise_key(key)
            child_path = path + "." + str(key) if path else str(key)
            if normalized in _COEFFICIENT_KEY_ALIASES or "coefficient" in normalized:
                found.append({
                    "path": child_path,
                    "key": str(key),
                    "shape": _shape(value),
                })
            if isinstance(value, Mapping):
                visit(value, child_path, depth + 1)

    visit(archive, "", 0)
    return found[:128]


def _structural_comparison(
    loaded: Any,
    archives: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    fields = ("M", "M_Z", "M_inverse", "y", "y_Z")
    actual_shapes = {name: _shape(getattr(loaded, name)) for name in fields}
    actual_coefficients = getattr(loaded, "coefficients", None)
    actual_coefficient_keys = (
        sorted(str(key) for key in actual_coefficients.keys())
        if isinstance(actual_coefficients, Mapping)
        else None
    )
    graph = getattr(loaded, "graph")
    graph_key, graph_rows = _graph_rows(graph)
    graph_summary = _rows_summary(graph_rows) if graph_rows is not None else None

    archive_reports: dict[str, Any] = {}
    for role, archive in archives.items():
        archive_key, archive_rows = _first_archive_rows(archive)
        archive_summary = _rows_summary(archive_rows) if archive_rows is not None else None
        archive_reports[role] = {
            "archive_graph_rows_key": archive_key,
            "archive_graph_rows": archive_summary,
            "coefficient_candidates": _coefficient_candidates(archive),
        }

    comparisons: list[dict[str, Any]] = []
    aliases = {
        "m": "M",
        "mz": "M_Z",
        "minverse": "M_inverse",
        "matrix": "M",
        "matrixz": "M_Z",
        "matrixinverse": "M_inverse",
        "y": "y",
        "yz": "y_Z",
    }
    for role, report in archive_reports.items():
        for candidate in report["coefficient_candidates"]:
            field = aliases.get(_normalise_key(candidate["key"]))
            if field is None:
                continue
            comparisons.append({
                "role": role,
                "archive_path": candidate["path"],
                "actual_field": field,
                "shape_equal": _shape_equal(candidate["shape"], actual_shapes[field]),
            })

    archive_coefficient_keys: dict[str, list[str] | None] = {}
    for role, archive in archives.items():
        bindings = archive.get("full_defect_linear_bindings")
        coefficients = bindings.get("coefficients") if isinstance(bindings, Mapping) else None
        archive_coefficient_keys[role] = (
            sorted(str(key) for key in coefficients.keys())
            if isinstance(coefficients, Mapping)
            else None
        )

    archive_counts = {
        role: report["archive_graph_rows"]["count"]
        if report["archive_graph_rows"] is not None else None
        for role, report in archive_reports.items()
    }
    graph_count = graph_summary["count"] if graph_summary is not None else None
    archive_boundaries = {
        role: {
            "first_row": report["archive_graph_rows"].get("first_row")
            if report["archive_graph_rows"] is not None else None,
            "last_row": report["archive_graph_rows"].get("last_row")
            if report["archive_graph_rows"] is not None else None,
            "operation_prefix": report["archive_graph_rows"].get("operation_prefix")
            if report["archive_graph_rows"] is not None else None,
            "operation_suffix": report["archive_graph_rows"].get("operation_suffix")
            if report["archive_graph_rows"] is not None else None,
        }
        for role, report in archive_reports.items()
    }
    graph_boundaries = {
        "first_row": graph_summary.get("first_row") if graph_summary is not None else None,
        "last_row": graph_summary.get("last_row") if graph_summary is not None else None,
        "operation_prefix": graph_summary.get("operation_prefix") if graph_summary is not None else None,
        "operation_suffix": graph_summary.get("operation_suffix") if graph_summary is not None else None,
    }
    return {
        "actual_fields": actual_shapes,
        "actual_graph": {
            "type": type(graph).__module__ + "." + type(graph).__qualname__,
            "rows_attribute": graph_key,
            "rows": graph_summary,
        },
        "archives": archive_reports,
        "coefficient_shape_comparisons": comparisons,
        "coefficient_key_comparison": {
            "actual_keys": actual_coefficient_keys,
            "archive_keys": archive_coefficient_keys,
            "reserved_keys_equal": (
                actual_coefficient_keys is not None
                and archive_coefficient_keys.get("reserved") == actual_coefficient_keys
            ),
            "full_keys_equal": (
                actual_coefficient_keys is not None
                and archive_coefficient_keys.get("full") == actual_coefficient_keys
            ),
        },
        "graph_row_count_comparison": {
            "actual_graph_count": graph_count,
            "archive_counts": archive_counts,
            "reserved_count_equal": graph_count is not None and graph_count == archive_counts.get("reserved"),
            "full_count_equal": graph_count is not None and graph_count == archive_counts.get("full"),
        },
        "graph_row_structure_comparison": {
            "actual_graph_boundaries": graph_boundaries,
            "archive_boundaries": archive_boundaries,
            "bounded_only": True,
            "complete_graph_scan": False,
            "reason": (
                "the restored consumer already verifies the exact reserved prefix "
                "and archive node relationships; this runner records no redundant full scan"
            ),
        },
        "structural_only": True,
        "source_values_evaluated": False,
    }


def _loaded_action_output(loaded: Any, progress: dict[str, Any]) -> dict[str, Any]:
    """Record the one action performed inside ``restore_actual_graph``."""

    action = getattr(loaded, "action", None)
    _require(action is not None, "restored graph did not retain its MatrixC1Action")
    progress["inverse_c1_action_calls"] = 1
    progress["direct_inverse_c1_action_calls"] = 0
    return {
        "action_call_owner": "restore_actual_graph",
        "action_call_count": 1,
        "inverse_action_report": _report_record(
            action.report,
            (
                "schema",
                "dimension",
                "same_graph_handles",
                "source_handles_bound",
                "caller_inverse_assumed",
                "inverse_identity_verified",
                "actual_source_owner_authenticated",
                "status",
                "derivative_formula",
            ),
        ),
        "beta_handles": _opaque_record(action.beta),
        "beta_Z_handles": _opaque_record(action.beta_Z),
        "beta_pair_handles": _opaque_record(action.beta_pairs),
    }


def _consume(
    hashes: Mapping[str, str],
    archives: Mapping[str, tuple[str, str, str]],
    progress: dict[str, Any],
) -> dict[str, Any]:
    archive_raw = {
        role: _recheck_raw(hashes, declaration[0])
        for role, declaration in archives.items()
    }
    progress["raw_recheck_complete"] = True

    # The graph consumer uses a package-relative import for the injected C1
    # action.  Import the pinned action under its package name first, then the
    # graph consumer under the same package namespace; neither module is
    # imported before the exclusive entry receipt.
    action_module = _import_new(C1_ACTION, hashes)
    progress["c1_action_consumer_imported"] = True
    actual_module = _import_new(ACTUAL_GRAPH, hashes)
    progress["actual_graph_consumer_imported"] = True
    _require(getattr(action_module, "SCHEMA", None) == "reserved_matrix_c1_action_v1",
             "injected C1 action schema changed")
    loader = getattr(actual_module, "restore_actual_graph", None)
    if not callable(loader):
        loader = getattr(actual_module, "load_reserved_actual_graph", None)
    _require(callable(loader), "actual graph loader is missing")

    # Decode only after raw recheck and entry.  The new consumer accepts
    # already-decoded mappings plus the externally verified raw SHA pins.
    decoded = {
        role: _load_json(gzip.decompress(raw), role + " archive")
        for role, raw in archive_raw.items()
    }
    _require(all(isinstance(value, Mapping) for value in decoded.values()),
             "decoded archive must be a JSON object")
    progress["gzip_decoded_after_entry"] = True

    loaded = loader(
        reserved_archive=decoded["reserved"],
        full_archive=decoded["full"],
        reserved_raw_sha256=hashes[archives["reserved"][0]],
        full_raw_sha256=hashes[archives["full"][0]],
        raw_sha256_verified=True,
        action=action_module.inverse_c1_action,
    )
    progress["actual_graph_loader_calls"] = 1
    for field in ("graph", "M", "M_Z", "M_inverse", "y", "y_Z"):
        _require(hasattr(loaded, field), "actual graph result lacks ." + field)

    handles = {
        field: _opaque_record(getattr(loaded, field))
        for field in ("M", "M_Z", "M_inverse", "y", "y_Z")
    }
    source_id_constants = {}
    for name in (
        "M_NODE_IDS",
        "INVERSE_NODE_IDS",
        "Y_NODE_IDS",
        "Y_Z_NODE_IDS",
        "COEFFICIENT_NODE_IDS",
    ):
        value = getattr(actual_module, name, _MISSING)
        if value is not _MISSING:
            source_id_constants[name] = _source_id_literal(value)
    structural = _structural_comparison(loaded, decoded)
    action = _loaded_action_output(loaded, progress)
    embedded_match = getattr(loaded, "coefficient_match", None)
    embedded_match_record: dict[str, Any] | None = None
    if embedded_match is not None:
        checked = getattr(embedded_match, "checked", None)
        methods = getattr(embedded_match, "methods", None)
        embedded_match_record = {
            "checked": dict(checked) if isinstance(checked, Mapping) else None,
            "methods": dict(methods) if isinstance(methods, Mapping) else None,
            "complete": getattr(embedded_match, "complete", None),
            "status": getattr(embedded_match, "status", None),
        }
    coefficient_handles = None
    if isinstance(getattr(loaded, "coefficients", None), Mapping):
        coefficient_handles = {
            str(key): _opaque_record(value)
            for key, value in getattr(loaded, "coefficients").items()
        }
    dependency_proofs = {}
    for name in ("mu_dependency", "matrix_dependency", "inverse_dependency"):
        proof = getattr(loaded, name, None)
        if proof is not None:
            dependency_proofs[name] = _report_record(
                proof, ("roots", "node_count", "operations", "Z_independent")
            )
    return {
        "actual_graph_handles": handles,
        "source_id_constants": source_id_constants,
        "actual_coefficient_handles": coefficient_handles,
        "embedded_coefficient_match": embedded_match_record,
        "dependency_proofs": dependency_proofs,
        "restored_graph_flags": {
            "reserved_prefix_verified": getattr(loaded, "reserved_prefix_verified", None),
            "source_owner": getattr(loaded, "source_owner", None),
            "actual_source_owner_authenticated": getattr(loaded, "actual_source_owner_authenticated", None),
            "inverse_identity_verified": getattr(loaded, "inverse_identity_verified", None),
        },
        "structural_comparison": structural,
        "c1_action": action,
        "archive_raw_sha256": {
            role: hashes[declaration[0]] for role, declaration in archives.items()
        },
        "claims": {
            "physical_closure": False,
            "inverse_identity_verified": False,
            "source_values_evaluated": False,
            "old_producer_imported": False,
            "fixture_called": False,
            "historical_ancestor_graph_traversed": False,
            "new_graph_dependency_closure_traversed": True,
            "numerical_evaluation": False,
        },
    }


def _error_record(exc: BaseException) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def main() -> int:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    if ENTRY.exists() or RESULT.exists():
        try:
            raise RunnerInputError("component already entered")
        except BaseException as exc:
            error = _error_record(exc)
        print(json.dumps({"component": NAME, "status": "rejected", "error": error}, sort_keys=True))
        return 1

    source_sha: str | None = None
    hashes: dict[str, str] = {}
    archives: dict[str, tuple[str, str, str]] = {}
    preflight_error: dict[str, str] | None = None
    try:
        _, hashes, archives, source_sha = _preflight()
    except BaseException as exc:
        preflight_error = _error_record(exc)

    entry_value: dict[str, Any] = {
        "component": NAME,
        "entered": True,
        "source_sha256": source_sha,
        "input_sha256": dict(sorted(hashes.items())),
        "archive_roles": {
            role: {"path": value[0], "raw_sha256": value[1]}
            for role, value in sorted(archives.items())
        },
        "preflight_passed": preflight_error is None,
        "raw_sha256_verified_before_entry": preflight_error is None,
        "gzip_decoded_before_entry": False,
        "old_producer_imported": False,
        "fixture_called": False,
        "numerical_evaluation": False,
    }
    if preflight_error is not None:
        entry_value["error"] = preflight_error
    try:
        entry_sha = _write_exclusive(ENTRY, entry_value)
    except BaseException as exc:
        print(json.dumps({"component": NAME, "status": "failed", "error": _error_record(exc)}, sort_keys=True))
        return 1

    progress: dict[str, Any] = {
        "raw_recheck_complete": False,
        "actual_graph_consumer_imported": False,
        "c1_action_consumer_imported": False,
        "gzip_decoded_after_entry": False,
        "actual_graph_loader_calls": 0,
        "inverse_c1_action_calls": 0,
        "direct_inverse_c1_action_calls": 0,
        "full_third_row_c1_action_calls": 0,
    }
    outcome: dict[str, Any] = {
        "component": NAME,
        "status": "failed",
        "entry_sha256": entry_sha,
        "source_sha256": source_sha,
        "input_sha256": dict(sorted(hashes.items())),
        "archive_roles": {
            role: {"path": value[0], "raw_sha256": value[1]}
            for role, value in sorted(archives.items())
        },
        "raw_sha256_verified_before_entry": preflight_error is None,
        "old_producer_imported": False,
        "fixture_called": False,
        "numerical_evaluation": False,
        "physical_closure_claimed": False,
        "progress": progress,
        "output": None,
    }
    if preflight_error is not None:
        outcome["error"] = preflight_error
    else:
        try:
            outcome["output"] = _consume(hashes, archives, progress)
            outcome["status"] = "returned"
        except BaseException as exc:
            outcome["error"] = _error_record(exc)
        finally:
            outcome["progress"] = progress

    try:
        result_sha = _write_exclusive(RESULT, outcome)
    except BaseException as exc:
        print(json.dumps({"component": NAME, "status": "failed", "entry_sha256": entry_sha,
                          "error": _error_record(exc)}, sort_keys=True))
        return 1
    print(json.dumps({
        "component": NAME,
        "status": outcome["status"],
        "entry_sha256": entry_sha,
        "result_sha256": result_sha,
        "error": (outcome.get("error") or {}).get("message"),
    }, sort_keys=True))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
