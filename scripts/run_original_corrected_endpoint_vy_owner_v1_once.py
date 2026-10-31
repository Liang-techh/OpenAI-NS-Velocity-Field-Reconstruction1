"""One-entry frozen consumer for the corrected endpoint ``V_y`` source owner.

All declared raw bytes are authenticated before the exclusive entry receipt is
written and again before the endpoint archive is decoded.  The synthetic
fixture is run independently after entry and before the real endpoint is
read.  The real consumer only retains the four signed source cells and their
formal scales; it does not run the endpoint producer, the corrected m15
producer, a graph walk, a numerical evaluator, or a receipt ancestor.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "original_corrected_endpoint_vy_owner_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

ENDPOINT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1.json.gz"
)
ENDPOINT_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1_check.json"
)
CORRECTED_RESERVED = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json"
)
ENDPOINT_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1.py"
)
CORRECTED_RESERVED_RUNNER = "scripts/run_original_corrected_reserved_m15.py"

INPUT = (
    ENDPOINT,
    ENDPOINT_CHECK,
    CORRECTED_RESERVED,
    ENDPOINT_SOURCE,
    CORRECTED_RESERVED_RUNNER,
)
LEAF = "src/openai_ns_reconstruction/original_corrected_endpoint_vy_owner_v1.py"
FIXTURE = "src/openai_ns_reconstruction/original_corrected_endpoint_vy_owner_fixture_v1.py"
RUNNER = "scripts/run_original_corrected_endpoint_vy_owner_v1_once.py"
REQUIRED = INPUT + (LEAF, FIXTURE, RUNNER)

FIXED = {
    ENDPOINT: "e691db8090b4400537228db14d387a46e3163e73c94880a319a04f576cf8ec78",
    ENDPOINT_CHECK: "02dec83e6e5b4c52d8aadc7e5f8dc359b2322dc62e5aa3ce895f60e8f855ecb2",
    CORRECTED_RESERVED: "fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5",
    ENDPOINT_SOURCE: "828afe7bd79f91d8ccbeb989711147306b41abe400344b4b9fa3f55aaa14eb37",
    CORRECTED_RESERVED_RUNNER: "9b41d60ab277656734cc98aa3272123a1ce383da4337ba85982dd41142350d6d",
}

CELLS = (("-1", "-.5"), ("-.5", "0"), ("0", ".5"), (".5", "1"))
RECEIPT_CHECKS = (
    "exact_four_Z_cell_cover",
    "lower_actual_V_y_rows_and_U_ZZ_checked",
    "physical_common_Pstar_once",
)


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
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        path = row.get("path")
        expected = row.get("sha256")
        require(isinstance(path, str) and path in REQUIRED, "unexpected manifest input")
        require(path not in hashes, "duplicate manifest input: " + path)
        require(isinstance(expected, str) and len(expected) == 64, "invalid manifest SHA: " + path)
        actual = sha(repository_path(path).read_bytes())
        require(actual == expected.lower(), "manifest SHA changed: " + path)
        if path in FIXED:
            require(actual == FIXED[path], "fixed SHA changed: " + path)
        hashes[path] = actual
    require(set(hashes) == set(REQUIRED), "manifest does not pin exactly the required inputs")
    for path in (ENDPOINT_SOURCE, LEAF, FIXTURE):
        require(path in hashes, "source file is not manifest-pinned: " + path)
    return hashes


def preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    require(BRIEF.is_file() and SOURCE.is_file(), "pending brief or source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = load_json(brief_raw, "pending brief")
    source_raw = SOURCE.read_bytes()
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
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "pending source checkpoint is missing")
    require(Path(str(checkpoint.get("path"))).resolve() == SOURCE.resolve(), "source checkpoint path changed")
    require(sha(source_raw) == checkpoint.get("sha256"), "source checkpoint SHA changed")
    require(source.get("component") == COMPONENT and source.get("one_entry_only") is True,
            "source manifest guard")
    hashes = manifest_hashes(source)
    return source, hashes, sha(source_raw), sha(brief_raw)


def authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    """Read and authenticate every declared raw byte before endpoint decode."""

    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = repository_path(path).read_bytes()
        require(sha(data) == hashes[path], "raw input changed after entry: " + path)
        raw[path] = data
    return raw


def module(relative: str, name: str) -> Any:
    path = repository_path(relative)
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "module load failed: " + relative)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def decode_endpoint(raw: bytes) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError("endpoint archive is not valid gzip") from exc
    return load_json(decoded, ENDPOINT + " decoded JSON")


def _endpoint_report(value: Mapping[str, Any]) -> Mapping[str, Any]:
    candidate = value.get("endpoint_report")
    if isinstance(candidate, Mapping):
        return candidate
    candidate = value.get("report")
    if isinstance(candidate, Mapping) and "candidate_preRc_endpoint_C1_constructed" in candidate:
        return candidate
    return value


def _receipt(value: Mapping[str, Any]) -> Mapping[str, Any]:
    if value.get("all_passed") is True:
        return value
    for key in ("endpoint_receipt", "receipt", "checker", "report_checker", "check"):
        candidate = value.get(key)
        if isinstance(candidate, Mapping):
            return candidate
    raise RunnerInputError("endpoint checker receipt is missing")


def _cell_cover(value: Any, label: str) -> tuple[tuple[str, str], ...]:
    require(isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
            label + " must be a cell sequence")
    result: list[tuple[str, str]] = []
    for index, item in enumerate(value):
        require(isinstance(item, Sequence) and not isinstance(item, (str, bytes, bytearray))
                and len(item) == 2 and all(isinstance(part, str) for part in item),
                f"{label}[{index}] is not an exact two-string cell")
        result.append((item[0], item[1]))
    require(tuple(result) == CELLS, label + " is not the canonical four-cell cover")
    return tuple(result)


def _source_row_flags(owner: Mapping[str, Any]) -> dict[str, Any]:
    cells = owner.get("cell_sources")
    require(isinstance(cells, Mapping), "owner cell_sources missing")
    retained_rows = 0
    selected_rows = 0
    for cell in CELLS:
        entry = cells.get(str(cell))
        require(isinstance(entry, Mapping), "owner cell missing: " + str(cell))
        lower = entry.get("lower")
        require(isinstance(lower, Mapping), "owner lower source missing: " + str(cell))
        source_rows = lower.get("source_rows")
        require(isinstance(source_rows, Mapping), "owner lower source_rows missing: " + str(cell))
        rows = source_rows.get("V_y")
        selected = source_rows.get("V_y_C0_C1")
        require(isinstance(rows, list) and len(rows) == 4, "four V_y rows were not retained")
        require(isinstance(selected, list) and len(selected) == 2, "V_y C0/C1 rows were not retained")
        for retained in rows:
            require(isinstance(retained, Mapping) and retained.get("raw_row_preserved") is True,
                    "raw signed source row was not preserved")
            require(retained.get("log_absolute_upper_not_used_as_value") is True
                    and retained.get("midpoint_not_used_as_value") is True,
                    "bound or midpoint was promoted to a source value")
            row = retained.get("row")
            require(isinstance(row, Mapping)
                    and row.get("encloses_original_source_function") is True
                    and row.get("point_value_selected") is False,
                    "V_y row source-enclosure flags changed")
            signed = retained.get("signed_coefficient")
            require(isinstance(signed, Mapping) and signed.get("source_key") is not None,
                    "signed coefficient was not retained")
            scaling = retained.get("symbolic_scaling")
            require(
                (isinstance(scaling, Mapping) and scaling.get("source_key") is not None)
                or any(
                    key in row
                    for key in ("symbolic_scaling", "formal_scaling", "source_scaling",
                                "scaling", "formal_positive_scale")
                ),
                "formal scaling was not retained",
            )
            retained_rows += 1
        selected_rows += len(selected)
    return {"retained_V_y_rows": retained_rows, "retained_C0_C1_rows": selected_rows}


def _corrected_scope(corrected: Mapping[str, Any], family: Any) -> Mapping[str, Any]:
    checks = corrected.get("checks")
    require(isinstance(checks, Mapping) and checks.get("passed") is True,
            "corrected reserved report checks.passed is not true")
    require(checks.get("old_212_six_reserved_counts_preserved") is True,
            "corrected reserved report changed the 212-term ledger")
    source_result = corrected.get("corrected_reserved_m15_source")
    require(isinstance(source_result, Mapping)
            and source_result.get("passed") is True
            and source_result.get("source_only") is True
            and source_result.get("source_executed") is False
            and source_result.get("existing_ledger_edited") is False
            and source_result.get("own_swirl_added") is False,
            "corrected reserved source scope changed")
    require(corrected.get("old_212_terms_and_all_errors_unchanged") is True
            and corrected.get("extra_reserved_term_added") is False,
            "corrected reserved ledger preservation changed")
    terminal = corrected.get("corrected_terminal_source_inputs")
    require(isinstance(terminal, Mapping) and terminal.get("source_only") is True,
            "corrected terminal source inputs are missing")
    require(terminal.get("numeric_endpoint_cell") == ["0", ".5"],
            "corrected terminal endpoint cell changed")
    proof = terminal.get("physical_endpoint_source_proof")
    require(isinstance(proof, Mapping) and proof.get("source_family") == family,
            "corrected terminal source family changed")
    return terminal


def consume(hashes: Mapping[str, str], stats: dict[str, Any]) -> dict[str, Any]:
    raw = authenticated_bytes(hashes)

    leaf = module(LEAF, "original_corrected_endpoint_vy_owner_v1")
    fixture = module(FIXTURE, "original_corrected_endpoint_vy_owner_fixture_v1")
    fixture_runner = getattr(fixture, "run_fixture", None)
    owner_builder = getattr(leaf, "build_owner", None)
    require(callable(fixture_runner) and callable(owner_builder), "leaf/fixture API is missing")

    fixture_result = fixture_runner(owner_builder)
    stats["fixture_invocations"] = 1
    require(isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True,
            "synthetic endpoint owner fixture failed")
    require(fixture_result.get("synthetic_only") is True
            and fixture_result.get("receipt_read") is False
            and fixture_result.get("report_file_read") is False
            and fixture_result.get("saved_source_graph_evaluated") is False
            and fixture_result.get("synthetic_formula_tree_evaluated") is True
            and fixture_result.get("midpoint_values_used") is False
            and fixture_result.get("log_absolute_upper_used_as_value") is False,
            "synthetic fixture scope changed")

    endpoint = decode_endpoint(raw[ENDPOINT])
    stats["archive_decode_started"] = True
    endpoint_check_outer = load_json(raw[ENDPOINT_CHECK], ENDPOINT_CHECK)
    corrected = load_json(raw[CORRECTED_RESERVED], CORRECTED_RESERVED)

    report = _endpoint_report(endpoint)
    receipt = _receipt(endpoint_check_outer)
    family = report.get("source_family")
    require(family is not None, "endpoint source_family is missing")
    require(receipt.get("all_passed") is True, "endpoint checker all_passed is not true")
    require(receipt.get("source_family") == family, "endpoint checker source_family changed")
    _cell_cover(report.get("cells"), "endpoint report.cells")
    _cell_cover(receipt.get("cells"), "endpoint receipt.cells")
    checks = receipt.get("checks")
    require(isinstance(checks, Mapping), "endpoint checker checks are missing")
    for name in RECEIPT_CHECKS:
        require(checks.get(name) is True, "endpoint checker missing " + name)

    terminal = _corrected_scope(corrected, family)
    samples = report.get("samples")
    require(isinstance(samples, Mapping), "endpoint report.samples is missing")
    lower = samples.get(str(("0", ".5")))
    require(isinstance(lower, Mapping) and isinstance(lower.get("lower"), Mapping),
            "endpoint report (0,.5) lower source is missing")
    saved_lower = lower["lower"]
    terminal_rows = terminal.get("physical_endpoint_source_rows")
    terminal_proof = terminal.get("physical_endpoint_source_proof")
    require(isinstance(terminal_rows, Mapping) and isinstance(terminal_proof, Mapping),
            "corrected terminal lower source rows or proof is missing")
    require(terminal_rows == saved_lower.get("source_rows"),
            "corrected terminal lower source rows differ from endpoint (0,.5)")
    require(terminal_proof == saved_lower.get("proof"),
            "corrected terminal lower source proof differs from endpoint (0,.5)")
    require(terminal_rows.get("V_y") == saved_lower["source_rows"].get("V_y"),
            "corrected (0,.5) same-source V_y rows are not identical")

    owner = owner_builder(
        endpoint,
        endpoint_receipt=receipt,
        expected_source_family=family,
        corrected_reserved_report=corrected,
        frozen=True,
        entered=True,
    )
    stats["owner_invocations"] = 1
    require(isinstance(owner, Mapping)
            and owner.get("schema") == "original_corrected_endpoint_vy_owner_v1"
            and owner.get("component") == "corrected_reserved_lower_endpoint_actual_physical_V_y"
            and owner.get("frozen_input_attested_by_caller") is True
            and owner.get("entered_input_attested_by_caller") is True,
            "endpoint owner returned an unexpected packet")
    require(owner.get("source_family") == family
            and owner.get("cells") == [list(cell) for cell in CELLS],
            "endpoint owner source family or cells changed")
    require(owner.get("ordinary_Z_orders") == [0, 1]
            and owner.get("ordinary_Z1_is_fixed_coordinate_derivative") is True
            and owner.get("ordinary_Z1_is_not_a_second_factorial_conversion") is True,
            "ordinary Z contract changed")
    row_counts = _source_row_flags(owner)
    comparison = owner.get("corrected_terminal_source_comparison")
    require(isinstance(comparison, Mapping)
            and comparison.get("checked") is True
            and comparison.get("physical_endpoint_source_rows_match") is True
            and comparison.get("physical_endpoint_source_proof_match") is True,
            "corrected terminal source comparison was not admitted")
    require(owner.get("numeric_source_values_selected") is False
            and owner.get("graph_walked") is False
            and owner.get("ancestor_imported") is False
            and owner.get("receipt_decoded_here") is False,
            "owner scope changed to numerical or ancestor evaluation")
    limits = owner.get("limits")
    require(isinstance(limits, Mapping)
            and limits.get("log_absolute_upper_as_function_value") is False
            and limits.get("reported_sign_or_midpoint_as_function_value") is False
            and limits.get("old_error_removed") is False,
            "owner promoted a bound or removed an error")

    return {
        "schema": COMPONENT + "_output_v1",
        "source_family": family,
        "cells": [list(cell) for cell in CELLS],
        "fixture": fixture_result,
        "owner": owner,
        "source_row_counts": row_counts,
        "same_source_V_y_0_.5": True,
        "signed_coefficients_and_formal_scales_retained": True,
        "log_caps_not_values": True,
        "m15_added": False,
        "m15_computed": False,
        "ledger_edited": False,
        "old_errors_removed": False,
        "global_closure": False,
        "all_four_cells_checked": True,
        "actual_endpoint_function_values_selected": False,
        "endpoint_archive_only_decoded_after_entry": True,
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
            "archive_decode_started": False,
            "fixture_invocations": 0,
            "owner_invocations": 0,
            "m15_added": False,
            "ledger_edited": False,
        },
    )
    stats = {
        "archive_decode_started": False,
        "fixture_invocations": 0,
        "owner_invocations": 0,
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
        "result": None,
        "error": None,
        "source_manifest": source,
    }
    try:
        outcome["result"] = consume(hashes, stats)
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
