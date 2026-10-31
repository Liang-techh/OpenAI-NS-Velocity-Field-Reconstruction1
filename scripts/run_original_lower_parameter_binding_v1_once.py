"""Run the frozen lower-parameter binding leaf exactly once.

The runner authenticates the pending brief and every manifest byte before it
writes an exclusive entry receipt.  The synthetic fixture is then run against
the leaf, without reading any real report.  Only after that fixture succeeds is
the previous four-cell receipt decoded and passed to the source-only leaf,
together with the native report, corrected report, and producer source text.

This consumer records exact saved ``delta`` and ``logPstar`` provenance only.
It does not claim a live endpoint identity, interval equality as function
equality, whole-Z closure, or removal of any historical error.
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
COMPONENT = "original_lower_parameter_binding_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

FOUR_CELL_RECEIPT = (
    "experiments/root_st073/original_corrected_fourcell_signed_repair_v1_receipt.json.gz"
)
NATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
CORRECTED = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json"
)
PRODUCER = "scripts/run_original_corrected_reserved_m15.py"
LEAF = "src/openai_ns_reconstruction/original_lower_parameter_binding_v1.py"
FIXTURE = "src/openai_ns_reconstruction/original_lower_parameter_binding_fixture_v1.py"
RUNNER = "scripts/run_original_lower_parameter_binding_v1_once.py"

FIXED = {
    FOUR_CELL_RECEIPT: "5876bf565709eab3e280b3e49b68cff57c30e5c772cb22a0c6c9f8b3f3b65b78",
    NATIVE: "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a",
    CORRECTED: "fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5",
    PRODUCER: "9b41d60ab277656734cc98aa3272123a1ce383da4337ba85982dd41142350d6d",
}
REQUIRED = (FOUR_CELL_RECEIPT, NATIVE, CORRECTED, PRODUCER, LEAF, FIXTURE, RUNNER)
NEW_CODE = (LEAF, FIXTURE, RUNNER)


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
        require(isinstance(expected, str) and len(expected) == 64,
                "invalid manifest SHA: " + str(path))
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
    require(source.get("component") == COMPONENT and source.get("one_entry_only") is True,
            "source manifest guard")
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "pending brief source checkpoint is missing")
    require(Path(str(checkpoint.get("path"))).resolve() == SOURCE.resolve(),
            "pending brief points to another source manifest")
    require(checkpoint.get("sha256") == sha(source_raw),
            "pending brief source checkpoint SHA changed")
    source_checkpoint = source.get("source_checkpoint")
    if source_checkpoint is not None:
        require(isinstance(source_checkpoint, Mapping), "source checkpoint is not an object")
        require(Path(str(source_checkpoint.get("path"))).resolve() == SOURCE.resolve()
                and source_checkpoint.get("sha256") == sha(source_raw),
                "source manifest checkpoint changed")
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
    """Inspect only top-level API names; never evaluate source text here."""

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
        raise RunnerInputError("four-cell receipt is not valid gzip") from exc
    return load_json(decoded, FOUR_CELL_RECEIPT + " decoded JSON")


def validate_fixture(result: Mapping[str, Any]) -> None:
    require(result.get("passed") is True, "synthetic parameter fixture failed")
    require(result.get("synthetic_only") is True
            and result.get("real_sources_read") is False
            and result.get("implementation_imported") is False,
            "synthetic fixture scope changed")
    checks = result.get("checks")
    require(isinstance(checks, Mapping), "synthetic fixture checks are missing")
    expected = (
        "accepted_exact_records",
        "direct_ast_producers",
        "live_identity_not_claimed",
        "z1_not_rederived",
        "rejects_native_hash_mismatch",
        "rejects_one_bit_endpoint_change",
        "rejects_incorrect_direct_producer",
    )
    require(all(checks.get(name) is True for name in expected),
            "synthetic fixture did not pass every guard")


def validate_leaf(result: Mapping[str, Any]) -> dict[str, Any]:
    require(result.get("schema") == "original_lower_parameter_binding_v1",
            "parameter leaf schema changed")
    require(result.get("passed") is True, "parameter leaf did not pass")
    scope = result.get("provenance_scope")
    require(isinstance(scope, Mapping), "parameter provenance scope is missing")
    for key in (
        "endpoint_live_owner_identity",
        "delta_min_branch_identity",
        "Pstar_source_function_identity",
        "cell_basis_identity",
        "native_corrected_endpoint_identity",
        "interval_equality_promoted_to_function_identity",
    ):
        require(scope.get(key) is False, "parameter leaf overclaimed " + key)
    require(result.get("source_execution") is False
            and result.get("real_source_imported") is False
            and result.get("receipt_decoded_here") is False
            and result.get("graph_traversed") is False,
            "parameter leaf execution scope changed")
    delta = result.get("delta")
    logpstar = result.get("actual_logPstar_parameter")
    require(isinstance(delta, Mapping) and delta.get("exact_saved_record_match") is True,
            "exact delta record match is missing")
    require(isinstance(logpstar, Mapping)
            and logpstar.get("exact_saved_record_match") is True,
            "exact logPstar record match is missing")
    repair = result.get("repair_receipt_contract")
    require(isinstance(repair, Mapping)
            and repair.get("source_family_matches") is True
            and repair.get("ledger_written") is False
            and repair.get("old_errors_removed") is False
            and repair.get("source_rows_are_not_point_values") is True
            and repair.get("endpoint_function_values_selected") is False,
            "four-cell repair contract scope changed")
    ast_contract = result.get("source_ast_contract")
    require(isinstance(ast_contract, Mapping)
            and ast_contract.get("parsed") is True
            and ast_contract.get("native_read_pinned_assignment") is True
            and ast_contract.get("delta_direct_deepcopy") is True
            and ast_contract.get("actual_logPstar_parameter_direct_deepcopy") is True,
            "producer source AST contract is incomplete")
    ledger = result.get("ledger_scope")
    require(isinstance(ledger, Mapping)
            and ledger.get("old_graphs_unchanged") is True
            and ledger.get("old_212_terms_and_all_errors_unchanged") is True
            and ledger.get("existing_ledger_edited") is False,
            "ledger preservation scope changed")
    ordinary = result.get("ordinary_Z1_scope")
    require(isinstance(ordinary, Mapping)
            and ordinary.get("source_rows_retained") is True
            and ordinary.get("second_derivative_or_second_factorial_conversion") is False
            and ordinary.get("product_rule_rederived_here") is False,
            "ordinary Z1 scope changed")
    return {
        "schema": result["schema"],
        "passed": True,
        "exact_delta_record_match": True,
        "exact_logPstar_record_match": True,
        "live_function_identity": False,
        "interval_equality_as_function_identity": False,
        "old_errors_removed": False,
        "whole_Z_closure": False,
    }


def consume(source: Mapping[str, Any], hashes: Mapping[str, str], raw: Mapping[str, bytes], stats: dict[str, Any]) -> dict[str, Any]:
    leaf_static = static_api(raw, LEAF, "bind")
    fixture_static = static_api(raw, FIXTURE, "run_fixture")
    leaf_module = module(LEAF, "original_lower_parameter_binding_v1")
    fixture_module = module(FIXTURE, "original_lower_parameter_binding_fixture_v1")
    bind = getattr(leaf_module, "bind", None)
    run_fixture = getattr(fixture_module, "run_fixture", None)
    require(callable(bind) and callable(run_fixture), "parameter leaf or fixture API missing")

    # The fixture receives only synthetic mappings and is run before real data is read.
    fixture_result = run_fixture(bind)
    stats["fixture_invocations"] = 1
    require(isinstance(fixture_result, Mapping), "fixture returned a non-object")
    validate_fixture(fixture_result)

    # Raw SHA authentication happened before entry and again after entry.
    # Only now is the previous compressed receipt decoded.
    stats["receipt_decode_started"] = True
    repair_receipt = decode_receipt(raw[FOUR_CELL_RECEIPT])
    native = load_json(raw[NATIVE], NATIVE)
    corrected = load_json(raw[CORRECTED], CORRECTED)
    producer_text = raw[PRODUCER].decode("utf-8")
    stats["real_source_reads_started"] = True

    stats["leaf_invocations"] = 1
    leaf_result = bind(native, corrected, repair_receipt, producer_text)
    require(isinstance(leaf_result, Mapping), "parameter leaf returned a non-object")
    leaf_summary = validate_leaf(leaf_result)
    return {
        "schema": COMPONENT + "_output_v1",
        "api_static": {"leaf": leaf_static, "fixture": fixture_static},
        "fixture": dict(fixture_result),
        "leaf": dict(leaf_result),
        "leaf_summary": leaf_summary,
        "source_provenance": {
            "native": {"path": NATIVE, "sha256": hashes[NATIVE]},
            "corrected": {"path": CORRECTED, "sha256": hashes[CORRECTED]},
            "producer": {"path": PRODUCER, "sha256": hashes[PRODUCER]},
            "four_cell_receipt": {"path": FOUR_CELL_RECEIPT, "sha256": hashes[FOUR_CELL_RECEIPT]},
            "leaf": {"path": LEAF, "sha256": hashes[LEAF]},
            "fixture": {"path": FIXTURE, "sha256": hashes[FIXTURE]},
        },
        "raw_sha_verified_before_receipt_decode": True,
        "cross_run_object_identity_assumed": False,
        "interval_equality_promoted_to_function_identity": False,
        "old_errors_removed": False,
        "whole_Z_closure": False,
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
        "entry_sha256": entry_sha,
        "result_sha256": result_sha,
    }, sort_keys=True))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
