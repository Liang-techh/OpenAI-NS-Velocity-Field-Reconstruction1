"""Run the frozen core-parameter lineage audit exactly once.

All manifest bytes are authenticated before the exclusive entry is written.
After entry, the new leaf and its synthetic fixture are loaded; the fixture
must pass before the previous compressed parameter receipt is decoded or the
native JSON is interpreted.  The leaf receives the decoded previous receipt
so its sibling ``source_provenance`` remains available for authentication.

The leaf's complete comparison result is preserved even when its source
status is ``mismatch`` or ``incomplete``.  Those statuses are evidence, not
runner failures and are never relabeled as proof.  This runner does not
import or execute either producer source, does not promote interval overlap
to function identity, and makes no live endpoint, basis, or native/corrected
RQy claim.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "original_core_parameter_lineage_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

NATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
PARAMETER_RECEIPT = (
    "experiments/root_st073/original_lower_parameter_binding_v1_receipt.json.gz"
)
Q_SOURCE = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
BACKGROUND_SOURCE = "src/openai_ns_reconstruction/original_micro_background_I1_numeric.py"
LEAF = "src/openai_ns_reconstruction/original_core_parameter_lineage_v1.py"
FIXTURE = "src/openai_ns_reconstruction/original_core_parameter_lineage_fixture_v1.py"
RUNNER = "scripts/run_original_core_parameter_lineage_v1_once.py"

PRODUCER_SOURCES = (Q_SOURCE, BACKGROUND_SOURCE)
REQUIRED = (NATIVE, PARAMETER_RECEIPT, Q_SOURCE, BACKGROUND_SOURCE, LEAF, FIXTURE, RUNNER)
NEW_CODE = (LEAF, FIXTURE, RUNNER)

FIXED = {
    NATIVE: "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a",
    PARAMETER_RECEIPT: (
        "355b1bfdfa371374bff547aa4e75185386053838f1a252b6d98bced28d104556"
    ),
}


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


def require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON constant is not admitted: " + token)


def load_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite JSON") from exc
    require(isinstance(value, dict), label + " must be a JSON object")
    return value


def save(path: Path, value: object) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    try:
        with path.open("xb") as handle:
            handle.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry receipt already exists: " + str(path)) from exc
    return sha(raw)


def repository_path(relative: str) -> Path:
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts, "non-relative input path")
    resolved = (ROOT / path).resolve()
    require(resolved.is_relative_to(ROOT.resolve()), "input path escapes repository: " + relative)
    require(resolved.is_file(), "missing input: " + relative)
    return resolved


def manifest_hashes(manifest: Mapping[str, Any]) -> dict[str, str]:
    rows = manifest.get("files")
    require(isinstance(rows, list), "source manifest files list is missing")
    allowed = set(REQUIRED)
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        path = row.get("path")
        expected = row.get("sha256")
        require(isinstance(path, str) and path in allowed, "unexpected manifest input")
        require(path not in hashes, "duplicate manifest input: " + path)
        require(
            isinstance(expected, str) and len(expected) == 64,
            "invalid manifest SHA: " + str(path),
        )
        actual = sha(repository_path(path).read_bytes())
        require(actual == expected.lower(), "manifest SHA changed: " + path)
        if path in FIXED:
            require(actual == FIXED[path], "fixed SHA changed: " + path)
        hashes[path] = actual
    require(set(REQUIRED) == set(hashes), "manifest must pin exactly the frozen inputs")
    require(all(path in hashes for path in NEW_CODE), "new code files are not all frozen")
    return hashes


def preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    require(BRIEF.is_file(), "pending brief is missing")
    require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = load_json(brief_raw, "pending brief")
    source = load_json(source_raw, "source manifest")
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == PHASE
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "pending brief guard",
    )
    require(
        source.get("component") == COMPONENT and source.get("one_entry_only") is True,
        "source manifest guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "pending brief source checkpoint is missing")
    require(Path(str(checkpoint.get("path"))).resolve() == SOURCE.resolve(),
            "pending brief points to another source manifest")
    require(checkpoint.get("sha256") == sha(source_raw),
            "pending brief source checkpoint SHA changed")
    source_checkpoint = source.get("source_checkpoint")
    if source_checkpoint is not None:
        require(isinstance(source_checkpoint, Mapping), "source checkpoint is not an object")
        require(
            Path(str(source_checkpoint.get("path"))).resolve() == SOURCE.resolve()
            and source_checkpoint.get("sha256") == sha(source_raw),
            "source manifest checkpoint changed",
        )
    hashes = manifest_hashes(source)
    return source, hashes, sha(source_raw), sha(brief_raw)


def authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for path in hashes:
        data = repository_path(path).read_bytes()
        require(sha(data) == hashes[path], "raw input changed after entry: " + path)
        raw[path] = data
    return raw


def static_api(raw: Mapping[str, bytes], path: str, function_name: str) -> dict[str, Any]:
    """Inspect only the new module's top-level API, without evaluating it."""

    source_text = raw[path].decode("utf-8")
    try:
        tree = ast.parse(source_text, filename=path)
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise RunnerInputError("source API is not parseable: " + path) from exc
    functions = {
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    require(function_name in functions, "required API is missing: " + path + ":" + function_name)
    return {"path": path, "function": function_name, "static_only": True}


def module(relative: str, name: str) -> Any:
    path = repository_path(relative)
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "module load failed: " + relative)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def decode_receipt(raw: bytes) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError("parameter-binding receipt is not valid gzip") from exc
    return load_json(decoded, PARAMETER_RECEIPT + " decoded JSON")


def validate_fixture(result: Mapping[str, Any]) -> None:
    require(result.get("passed") is True, "synthetic core-lineage fixture failed")
    require(
        result.get("synthetic_only") is True
        and result.get("real_sources_read") is False
        and result.get("implementation_imported") is False,
        "synthetic fixture scope changed",
    )
    checks = result.get("checks")
    require(isinstance(checks, Mapping), "synthetic fixture checks are missing")
    expected = (
        "accepts_exact_binary_lineage",
        "status_proved",
        "source_family_checked",
        "producer_ast_paths_checked",
        "display_metadata_not_identity",
        "overlap_not_used",
        "rejects_wrong_delta_read_path",
        "rejects_wrong_actual_record_binding",
        "rejects_one_bit_binary_change",
        "rejects_cross_source_family",
        "rejects_source_hash_mismatch",
    )
    require(all(checks.get(name) is True for name in expected),
            "synthetic fixture did not pass every guard")


def validate_leaf(result: Mapping[str, Any]) -> dict[str, Any]:
    require(result.get("schema") == "original_core_parameter_lineage_v1",
            "core-lineage leaf schema changed")
    status = result.get("status")
    require(status in {"proved", "mismatch", "incomplete"},
            "core-lineage source status is invalid")
    require(result.get("passed") is (status == "proved"),
            "core-lineage passed/status relation changed")
    require(
        result.get("source_execution") is False
        and result.get("real_receipt_opened_here") is False
        and result.get("producer_imported") is False
        and result.get("graph_traversed") is False,
        "core-lineage leaf execution scope changed",
    )
    provenance = result.get("provenance_scope")
    require(isinstance(provenance, Mapping), "core-lineage provenance scope is missing")
    for key in (
        "cross_run_object_identity_claimed",
        "live_endpoint_basis_identity",
        "native_corrected_RQy_identity",
        "source_function_point_value_identity",
        "interval_overlap_promoted_to_identity",
        "rounding_difference_promoted_to_identity",
    ):
        require(provenance.get(key) is False, "core-lineage overclaimed " + key)
    comparisons = result.get("comparisons")
    require(isinstance(comparisons, Mapping), "core-lineage comparisons are missing")
    require(comparisons.get("overlap_used_as_identity") is False,
            "interval overlap was used as identity")
    hashes = result.get("source_hashes")
    require(isinstance(hashes, Mapping), "core-lineage source hash checks are missing")
    for key in ("Q_increment_producer", "background_scale_producer"):
        item = hashes.get(key)
        require(isinstance(item, Mapping) and item.get("matches") is True,
                "native implementation hash did not authenticate " + key)
    return {
        "schema": result["schema"],
        "status": status,
        "passed": status == "proved",
        "mismatch_count": comparisons.get("mismatch_count"),
        "overlap_used_as_identity": False,
        "live_endpoint_basis_identity": False,
        "native_corrected_RQy_identity": False,
        "old_errors_removed": False,
    }


def consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    raw: Mapping[str, bytes],
    stats: dict[str, Any],
) -> dict[str, Any]:
    leaf_static = static_api(raw, LEAF, "verify")
    fixture_static = static_api(raw, FIXTURE, "run_fixture")
    leaf_module = module(LEAF, "original_core_parameter_lineage_v1")
    fixture_module = module(FIXTURE, "original_core_parameter_lineage_fixture_v1")
    verify = getattr(leaf_module, "verify", None)
    run_fixture = getattr(fixture_module, "run_fixture", None)
    require(callable(verify) and callable(run_fixture),
            "core-lineage leaf or fixture API missing")

    # The fixture receives only synthetic mappings and runs before real input
    # decoding.  This is the only leaf call made before the saved receipt is
    # opened.
    fixture_result = run_fixture(verify)
    stats["fixture_invocations"] = 1
    require(isinstance(fixture_result, Mapping), "fixture returned a non-object")
    validate_fixture(fixture_result)

    # Decode the prior receipt only after the synthetic fixture has passed.
    stats["receipt_decode_started"] = True
    previous_receipt = decode_receipt(raw[PARAMETER_RECEIPT])
    native = load_json(raw[NATIVE], NATIVE)
    q_source_text = raw[Q_SOURCE].decode("utf-8")
    background_source_text = raw[BACKGROUND_SOURCE].decode("utf-8")
    stats["real_source_reads_started"] = True

    stats["leaf_invocations"] = 1
    leaf_result = verify(
        native,
        previous_receipt,
        q_source_text,
        background_source_text,
    )
    require(isinstance(leaf_result, Mapping), "core-lineage leaf returned a non-object")
    leaf_summary = validate_leaf(leaf_result)
    source_status = leaf_summary["status"]
    return {
        "schema": COMPONENT + "_output_v1",
        "api_static": {"leaf": leaf_static, "fixture": fixture_static},
        "fixture": dict(fixture_result),
        "leaf": dict(leaf_result),
        "leaf_summary": leaf_summary,
        "source_comparison_status": source_status,
        "source_comparison_status_preserved": source_status in {
            "proved", "mismatch", "incomplete"
        },
        "source_comparison_mismatch_is_not_identity_proof": source_status != "proved",
        "previous_parameter_receipt_sha256": hashes[PARAMETER_RECEIPT],
        "native_source_sha256": hashes[NATIVE],
        "producer_source_provenance": {
            path: {"path": path, "sha256": hashes[path]}
            for path in PRODUCER_SOURCES
        },
        "source_provenance": {
            path: {"path": path, "sha256": hashes[path]}
            for path in (NATIVE, PARAMETER_RECEIPT, *PRODUCER_SOURCES, LEAF, FIXTURE)
        },
        "raw_sha_verified_before_receipt_decode": True,
        "actual_live_endpoint_owner_identity": False,
        "native_corrected_RQy_identity": False,
        "interval_equality_promoted_to_function_identity": False,
        "old_errors_removed": False,
        "ledger_edited": False,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = preflight()
    entry_sha = save(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "manifest_input_count": len(hashes),
            "receipt_decode_started": False,
            "fixture_invocations": 0,
            "leaf_invocations": 0,
            "old_errors_removed": False,
        },
    )
    stats: dict[str, Any] = {
        "receipt_decode_started": False,
        "real_source_reads_started": False,
        "fixture_invocations": 0,
        "leaf_invocations": 0,
    }
    outcome: dict[str, Any] = {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": "failed",
        "source_checkpoint_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "source_manifest": source,
        "result": None,
        "error": None,
        "old_errors_removed": False,
    }
    try:
        raw = authenticated_bytes(hashes)
        outcome["result"] = consume(source, hashes, raw, stats)
        outcome["status"] = "returned"
        outcome["source_comparison_status"] = outcome["result"]["source_comparison_status"]
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    outcome.update(stats)
    result_sha = save(RESULT, outcome)
    print(json.dumps({
        "component": COMPONENT,
        "phase": PHASE,
        "status": outcome["status"],
        "source_comparison_status": outcome.get("source_comparison_status"),
        "entry_sha256": entry_sha,
        "result_sha256": result_sha,
    }, sort_keys=True))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
