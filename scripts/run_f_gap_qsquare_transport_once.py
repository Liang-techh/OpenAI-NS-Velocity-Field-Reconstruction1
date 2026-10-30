"""Run the bounded F gap ``q_squared`` transport leaf exactly once.

This runner is the persistence boundary for
``f_gap_qsquare_transport_v1``.  It authenticates the seven raw inputs named
by the pending manifest before writing the exclusive-entry receipt and imports
only ``f_gap_qsquare_transport.py`` after that receipt.  The leaf owns the
single whole-JSON F-archive decode and may inspect only the two declared F
rows (``279381`` and ``279477``) plus at most two non-protected product-child
rows.  The runner never opens a prior mathematical result, follows graph
edges, runs BFS, replays an ancestor, evaluates a source, evaluates an
integral, or selects a numeric value.

Structural match and structural-unresolved packets are both valid returned
outcomes.  The formal correction recorded here is intentionally conditional:
it requires same-source ownership, capture-free substitution, and the
appropriate differentiability/ordinary-Z and moving-endpoint hypotheses.
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
COMPONENT = "f_gap_qsquare_transport_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_gap_qsquare_transport_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_gap_qsquare_transport.py")
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
PRIOR_QROLES_RESULT = Path("tmp/f_pressure_path_integral_roles_v1/result.json")
PRIOR_GAP_RESULT = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
NONRADIAL_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields.py"
)
PULSE_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py"
)

REQUIRED_INPUTS = (
    LEAF.as_posix(),
    RUNNER.as_posix(),
    F_ARCHIVE.as_posix(),
    PRIOR_QROLES_RESULT.as_posix(),
    PRIOR_GAP_RESULT.as_posix(),
    NONRADIAL_SOURCE.as_posix(),
    PULSE_SOURCE.as_posix(),
)

# The three archived inputs have fixed published pins.  Code and runner pins
# are supplied by the frozen manifest and are authenticated from raw bytes.
F_ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
PRIOR_QROLES_RESULT_SHA256 = (
    "ecbaeda54bdfa82952e09ed1639498ea0bb057716ffc087ab078aa31673c780b"
)
PRIOR_GAP_RESULT_SHA256 = (
    "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2"
)
EXPECTED_FIXED_SHA256 = {
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
    PRIOR_QROLES_RESULT.as_posix(): PRIOR_QROLES_RESULT_SHA256,
    PRIOR_GAP_RESULT.as_posix(): PRIOR_GAP_RESULT_SHA256,
}

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

PROTECTED_NODE_IDS = frozenset(
    {
        0,
        1,
        17,
        33572,
        25220,
        263344,
        263348,
        263349,
        262918,
        205147,
        263346,
    }
)
TARGET_NODE_IDS = frozenset({33982})
REQUIRED_DIRECT_ROWS = frozenset({279381, 279477})
MAX_DIRECT_ROWS = 4

# These are semantic conclusions that this bounded leaf may not claim.  The
# names are checked recursively so a nested role packet cannot silently
# promote a structural observation to a source or numerical proof.
FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "protected_dereferenced",
        "protected_ref_dereferenced",
        "protected_frontier_dereferenced",
        "protected_indexed",
        "protected_root_indexed",
        "protected_nodes_indexed",
        "target_node_indexed",
        "target_nodes_indexed",
        "target_row_indexed",
        "target_node_dereferenced",
        "target_row_dereferenced",
        "radius_row_indexed",
        "saved_target_dereferenced",
        "source_values_computed",
        "actual_source_values",
        "source_value_selected",
        "source_values_selected",
        "numeric_values_computed",
        "numeric_use_admissible",
        "integrals_evaluated",
        "integral_evaluated",
        "quadrature_evaluated",
        "graph_traversal",
        "graph_edge_search",
        "source_graph_traversed",
        "old_graph_traversed",
        "bfs_used",
        "ancestor_replayed",
        "old_ancestor_replayed",
        "old_leaf_imported",
        "old_producer_called",
        "duplicate_graph_computation",
        "repeated_graph_computation",
        "source_identity_proved",
        "same_source_function_proved",
        "function_identity_proved",
        "actual_source_identity_proved",
        "capture_free",
        "capture_free_proved",
        "capture_free_substitution_proved",
        "ordinary_Z_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "transport_proved",
        "endpoint_identity_proved",
        "endpoint_transport_proved",
        "moving_endpoint_rule_proved",
        "closure_proved",
        "q_square_closure_proved",
        "five_moment_closure_proved",
        "error_removal",
        "old_error_removed",
        "errors_removed",
    }
)

SEMANTIC_FALSE_KEYS = frozenset(
    {
        "source_identity_proved",
        "same_source_function_proved",
        "function_identity_proved",
        "capture_free",
        "capture_free_proved",
        "capture_free_substitution_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "transport_proved",
        "closure_proved",
        "q_square_closure_proved",
        "source_values_computed",
        "actual_source_values",
        "numeric_values_computed",
        "numeric_use_admissible",
        "error_removal",
        "old_error_removed",
        "errors_removed",
    }
)

DIRECT_ROW_KEYS = frozenset(
    {
        "direct_graph_row_ids",
        "indexed_graph_row_ids",
        "declared_graph_row_ids",
        "direct_row_ids",
        "indexed_row_ids",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    """Hash raw bytes only; do not decode an authenticated input here."""

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
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=repr)
    if isinstance(value, Path):
        return value.as_posix()
    return repr(value)


def _write_once(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(_json_safe(value))
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, dict), label + " is not a JSON object")
    return value


def _repository_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " is missing or escapes repository: " + normalized,
    )
    return key, path


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
        for field in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256"):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend(
                    (raw_path, expected, f"source.{field}")
                    for raw_path, expected in table.items()
                )
    _require(entries, "source manifest contains no direct raw SHA entries")
    return entries


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate exactly seven raw paths before any leaf import."""

    hashes: dict[str, str] = {}
    for raw_path, expected, label in _manifest_entries(source):
        key, path = _repository_path(raw_path, label)
        _require(key not in hashes, "duplicate manifest path: " + key)
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        digest = expected.lower()
        actual = _sha256_file(path)
        _require(actual.lower() == digest, "raw SHA changed: " + key)
        fixed = EXPECTED_FIXED_SHA256.get(key)
        if fixed is not None:
            _require(digest == fixed, "fixed SHA pin changed: " + key)
        hashes[key] = actual

    required = set(REQUIRED_INPUTS)
    _require(
        set(hashes) == required,
        "manifest must freeze exactly seven inputs; missing or extra paths: "
        + ", ".join(sorted((required ^ set(hashes)))),
    )
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(f"{COMPONENT}_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    _require(isinstance(raw_path, str) and raw_path, "source checkpoint path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(), "pending brief checkpoint points elsewhere")
    expected = checkpoint.get("sha256")
    actual = _sha256_file(SOURCE)
    _require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and actual.lower() == expected.lower(),
        "source manifest checkpoint SHA changed",
    )
    return actual


def _guard_true(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(any(source[name] is True for name in present), "manifest guard is not true: " + label)


def _guard_false(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(all(source[name] is False for name in present), "manifest guard is not false: " + label)


def _walk_true_keys(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = path + "." + key_text
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                found.append(child)
            found.extend(_walk_true_keys(item, child))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_true_keys(item, path + "[" + str(index) + "]"))
    return found


def _walk_key_values(value: Any, names: set[str], path: str = "$") -> list[tuple[str, Any]]:
    found: list[tuple[str, Any]] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = path + "." + key_text
            if key_text in names:
                found.append((child, item))
            found.extend(_walk_key_values(item, names, child))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_key_values(item, names, path + "[" + str(index) + "]"))
    return found


def _validate_false_scope(value: Mapping[str, Any], label: str) -> None:
    forbidden = _walk_true_keys(value)
    _require(not forbidden, label + " claims forbidden true scope: " + ", ".join(forbidden))
    for path, marker in _walk_key_values(value, set(SEMANTIC_FALSE_KEYS)):
        _require(marker is False, label + " semantic flag is true: " + path)


def _require_fixture(value: Any) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    passed = value.get("passed") is True or value.get("checks_passed") is True
    if not passed:
        checks = value.get("checks")
        passed = isinstance(checks, Mapping) and bool(checks) and all(item is True for item in checks.values())
    _require(passed, "synthetic fixture did not pass")
    _validate_false_scope(value, "synthetic fixture")
    _json_bytes(_json_safe(value))
    return value


def _row_ids(value: Any) -> list[int] | None:
    if isinstance(value, (list, tuple, set, frozenset)):
        if all(isinstance(item, int) and not isinstance(item, bool) for item in value):
            return list(value)
        ids: list[int] = []
        for item in value:
            if not isinstance(item, Mapping):
                return None
            candidate = next(
                (item[name] for name in ("node", "node_id", "row", "row_id", "id") if name in item),
                None,
            )
            if not isinstance(candidate, int) or isinstance(candidate, bool):
                return None
            ids.append(candidate)
        return ids
    return None


def _direct_row_scope(value: Mapping[str, Any]) -> dict[str, Any]:
    candidates: list[tuple[str, list[int]]] = []
    for path, marker in _walk_key_values(value, set(DIRECT_ROW_KEYS)):
        ids = _row_ids(marker)
        if ids is not None:
            candidates.append((path, ids))
    _require(candidates, "q-square result has no direct-row evidence")
    first = candidates[0][1]
    _require(len(first) == len(set(first)), "direct-row evidence contains duplicates")
    _require(REQUIRED_DIRECT_ROWS.issubset(first), "direct-row evidence omits 279381 or 279477")
    _require(len(first) <= MAX_DIRECT_ROWS, "q-square result indexed more than four F rows")
    _require(not (set(first) & (PROTECTED_NODE_IDS | TARGET_NODE_IDS)), "protected or target row was indexed")
    extras = sorted(set(first) - REQUIRED_DIRECT_ROWS)
    _require(len(extras) <= 2, "more than two product-child rows were indexed")
    for path, ids in candidates[1:]:
        _require(ids == first, "inconsistent direct-row evidence: " + path)
    return {
        "required_rows": sorted(REQUIRED_DIRECT_ROWS),
        "direct_rows": first,
        "product_child_rows": extras,
        "direct_row_count": len(first),
        "evidence_paths": [path for path, _ in candidates],
    }


def _import_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new q-square transport leaf is missing")
    spec = importlib.util.spec_from_file_location("f_gap_qsquare_transport_v1_new_leaf", path)
    _require(spec is not None and spec.loader is not None, "new q-square transport leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new leaf callable is missing: " + name)
    return function


FORMAL_CANDIDATE_CORRECTION = {
    "status": "conditional_formal_only",
    "formal_only": True,
    "coefficient": "c = -1/4",
    "definitions": {
        "M": "M = R^2/mu",
        "F": "F = M*Q0^2",
        "E_Q": "E_Q = Q1 - (partial_Z Q0)_xi",
    },
    "expression": (
        "D_Z(c*I0) - c*I1 = "
        "-1/4*integral_l^u(M_Z*Q0^2) dxi "
        "+1/2*integral_l^u(M*Q0*E_Q) dxi "
        "-1/4*F_u*u_Z + 1/4*F_l*l_Z"
    ),
    "premises_satisfied": False,
    "assumptions_unproved": [
        "same source owner for F integrand, R, mu, logRp, endpoints and measure",
        "capture-free old-to-F substitution preserving binders and operation tree",
        "ordinary-Z/fixed-xi differentiability for M and Q0, with Q1 source-bound",
        "upper-bound owner u = 13 - 4*mu_F and valid u_Z = -4*mu_F_Z",
        "lower-bound owner l = 11 and valid l_Z = 0",
        "moving-endpoint chain rule and endpoint evaluation of F = M*Q0^2",
        "no source values, integrals, numerical bounds, closure or error removal",
    ],
}


def _validate_transport(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "q-square transport leaf did not return an object")
    _validate_false_scope(value, "q-square transport leaf")
    row_scope = _direct_row_scope(value)
    return {
        "structural_only": True,
        "direct_row_scope": row_scope,
        "source_identity_proved": False,
        "capture_free_substitution_proved": False,
        "ordinary_Z_transport_proved": False,
        "closure_proved": False,
        "actual_source_values": False,
        "numeric_use_admissible": False,
        "error_removal": False,
    }


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief = _load_json_object(BRIEF, "pending brief")
    source = _load_json_object(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "brief component changed")
    _require(brief.get("phase") == PHASE, "brief phase changed")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal flags changed",
    )
    source_sha = _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "manifest component changed")
    _require(source.get("one_entry_only") is True, "manifest one-entry guard missing")
    _guard_true(
        source,
        (
            "all_raw_inputs_authenticated_before_decompression",
            "all_raw_inputs_authenticated_before_decode",
            "raw_sha256_authenticated_before_entry",
        ),
        "raw authentication",
    )
    _guard_true(
        source,
        ("archive_whole_json_decode", "archive_whole_json_decoded", "whole_json_archive_decode"),
        "whole archive JSON decode",
    )
    _guard_true(
        source,
        ("bounded_graph_direct_reads", "bounded_direct_graph_reads", "direct_graph_reads_bounded"),
        "bounded direct reads",
    )
    _guard_true(
        source,
        ("protected_and_target_nodes_not_individually_indexed", "no_protected_or_target_indexing"),
        "protected/target indexing",
    )
    _guard_true(
        source,
        ("unresolved_matches_preserved", "structural_match_or_unresolved_success"),
        "unresolved structural outcomes",
    )
    _guard_false(
        source,
        ("graph_traversal", "graph_edge_search", "source_graph_traversed", "bfs_used"),
        "graph traversal",
    )
    _guard_true(source, ("no_ancestor_replay", "no_old_ancestor_replay", "no_ancestor_math_replay"), "ancestor replay")
    _guard_true(source, ("no_integrals", "no_integrals_evaluated", "no_quadrature"), "integral evaluation")
    _guard_true(source, ("no_numeric_evaluation", "no_source_values", "no_numeric_values"), "numeric evaluation")
    hashes = _manifest_hashes(source)
    _require(LEAF.as_posix() in hashes, "new q-square leaf is not pinned")
    return dict(source), hashes, source_sha


def _recordable(value: Any) -> Any:
    return _json_safe(value) if value is not None else None


def main() -> int:
    # Unauthenticated pending metadata must not consume the one-entry slot.
    source, hashes, source_sha = _preflight()
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "new_leaf_imported": False,
            "consumer_calls_started": 0,
            "verification_invocations": 0,
            "archive_whole_json_decoded": False,
            "direct_graph_rows_indexed": [],
            "graph_traversal": False,
            "actual_source_values": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
        },
    )

    steps: list[dict[str, Any]] = []
    fixture_result: Any = None
    transport_result: Any = None
    transport_metrics: dict[str, Any] | None = None
    try:
        module = _import_leaf()
        fixture = _callable(module, "verify_synthetic_fixture")
        steps.append({"step": "verify_synthetic_fixture", "invoked": True})
        fixture_result = fixture()
        steps[-1]["returned_packet"] = _recordable(fixture_result)
        _require_fixture(fixture_result)

        builder = _callable(module, "build_qsquare_transport")
        steps.append({"step": "build_qsquare_transport", "invoked": True})
        # This is the sole source consumer call.  It receives only the root,
        # frozen manifest and already authenticated raw SHA mapping.
        transport_result = builder(ROOT, source, hashes)
        steps[-1]["returned_packet"] = _recordable(transport_result)
        transport_metrics = _validate_transport(transport_result)

        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": [item["step"] for item in steps],
            "verification_invocations": len(steps),
            "archive_whole_json_decoded": True,
            "direct_graph_rows_indexed": transport_metrics["direct_row_scope"]["direct_rows"],
            "graph_traversal": False,
            "actual_source_values": False,
            "source_values_computed": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "source_identity_proved": False,
            "capture_free_substitution_proved": False,
            "ordinary_Z_transport_proved": False,
            "closure_proved": False,
            "error_removal": False,
            "formal_candidate_correction": FORMAL_CANDIDATE_CORRECTION,
            "role_metrics": _json_safe(transport_metrics),
            "result": {
                "synthetic_fixture": _json_safe(fixture_result),
                "qsquare_transport": _json_safe(transport_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": None,
            "runtime_model_metadata": None,
        }
        _json_bytes(output)
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": [item["step"] for item in steps],
            "verification_invocations": len(steps),
            "archive_whole_json_decoded": bool(transport_metrics),
            "direct_graph_rows_indexed": [],
            "graph_traversal": False,
            "actual_source_values": False,
            "source_values_computed": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "source_identity_proved": False,
            "capture_free_substitution_proved": False,
            "ordinary_Z_transport_proved": False,
            "closure_proved": False,
            "error_removal": False,
            "formal_candidate_correction": FORMAL_CANDIDATE_CORRECTION,
            "role_metrics": _recordable(transport_metrics),
            "result": {
                "synthetic_fixture": _recordable(fixture_result),
                "qsquare_transport": _recordable(transport_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
            "runtime_model_metadata": None,
        }

    # Persist the complete JSON-safe packet before printing any status.
    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": output["consumer_calls"],
                "verification_invocations": output["verification_invocations"],
                "archive_whole_json_decoded": output["archive_whole_json_decoded"],
                "direct_graph_rows_indexed": output["direct_graph_rows_indexed"],
                "graph_traversal": output["graph_traversal"],
                "actual_source_values": output["actual_source_values"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
