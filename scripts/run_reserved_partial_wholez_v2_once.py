"""Consume the frozen reserved partial whole-Z source packet exactly once.

All frozen bytes are authenticated before the fixture or any gzip member is
decoded.  The independent fixture runs first; the new leaf is then called
once with exact interval inputs and no ancestor producer/helper.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "core_jet_projection_v2"
NAME = "reserved_partial_wholez_v2"
PHASE = NAME + "_pending"
DIRECTORY = ROOT / "tmp" / NAME
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF = "src/openai_ns_reconstruction/reserved_partial_wholez_v2.py"
FIXTURE = "src/openai_ns_reconstruction/reserved_partial_wholez_fixture_v2.py"
RUNNER = "scripts/run_reserved_partial_wholez_v2_once.py"
REPAIRED = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_repaired_Rp_m15_numeric.json"
)
COEFF = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_coefficient_enclosures.json.gz"
)
COEFF_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_coefficient_enclosures_check.json"
)
FULL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "full_five_numeric_bridge.json.gz"
)
FULL_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "full_five_numeric_bridge_check.json"
)
LINEAR = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_linear_operator.json.gz"
)
LINEAR_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_linear_operator_check.json"
)
FAILURE = "experiments/root_st073/reserved_partial_wholez_v1_failure.json.gz"
JOINT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_"
    "joint_reserved_m15_numeric.json"
)

FILES = (
    LEAF,
    FIXTURE,
    RUNNER,
    REPAIRED,
    COEFF,
    COEFF_CHECK,
    FULL,
    FULL_CHECK,
    LINEAR,
    LINEAR_CHECK,
    JOINT,
    FAILURE,
)
ARCHIVES = frozenset((COEFF, FULL, LINEAR, FAILURE))
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
INTERVAL_DPS = 130
FIXED_SHA256 = {
    FAILURE: "a1dc1f5ea4c1e613ac86a2061aea9ecbe1ca8ade8cf9c1996247d8ac179ef2ca",
    REPAIRED: "bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871",
    COEFF: "c4f9f43c52431195338cc1b06d6b399585ee47abefcda826806f51fe27600392",
    COEFF_CHECK: "a005a81c70a14e883b6bb365773027339bef09896398883929788cd1096072e2",
    FULL: "3dfebf4f49f7dd73f63c8418e0b18ff3e30c93a5860b4684d84df0620855a842",
    FULL_CHECK: "e00980861fbcfaca3bad0c8c0237ad0e1afb2ecd1df3163e91c6a7f17601dcea",
    LINEAR: "98fd45a48d72656ec91f19a5564577baed0dbb43fe85ce8c842336856dc17a17",
    LINEAR_CHECK: "b8fe16b23af26a51a3601c5c0dfdc888591d1a6303c74a8ba5d53da7c7b7827c",
    JOINT: "a46413fd755aa96c385aa7b81e1919fe0d22838868b8b5a7d4348708dd16dad6",
}


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON constant is not admitted: " + token)


def _load_json_bytes(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not valid finite JSON") from exc
    _require(isinstance(value, dict), label + " is not a JSON object")
    return value


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


def _json_safe(value: Any, *, context: Any = None, origin: str = "value") -> Any:
    if context is not None and hasattr(value, "_mpi_"):
        return _interval_packet(context, value, origin)
    if hasattr(value, "_mpf_"):
        return {
            "kind": "exact_mpf",
            "exact_mpf_tuple": list(tuple(value._mpf_)),
            "source_value_used": False,
        }
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {
            str(key): _json_safe(item, context=context, origin=origin)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_json_safe(item, context=context, origin=origin) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted(
            (_json_safe(item, context=context, origin=origin) for item in value),
            key=repr,
        )
    raise RunnerInputError("unsupported result value: " + type(value).__name__)


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(_json_safe(value))
    try:
        with path.open("xb") as stream:
            stream.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry receipt already exists") from exc
    return _sha256_bytes(raw)


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


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    rows = source.get("files")
    _require(isinstance(rows, list), "source manifest files missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        key, path = _repository_path(row.get("path"), f"source.files[{index}]")
        _require(key in FILES and key not in hashes, "unexpected or duplicate input: " + key)
        expected = row.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "manifest hash mismatch: " + key)
        fixed = FIXED_SHA256.get(key)
        if fixed is not None:
            _require(actual.lower() == fixed, "fixed source SHA changed: " + key)
        if key in ARCHIVES:
            _require(
                row.get("raw_sha256_before_decode") is True,
                "gzip input was not marked raw-authenticated: " + key,
            )
        hashes[key] = actual
    _require(set(hashes) == set(FILES), "manifest input set mismatch")
    return hashes


def _resolve_external_path(raw_path: Any) -> Path:
    _require(isinstance(raw_path, str) and raw_path, "checkpoint path missing")
    candidate = Path(raw_path)
    return candidate if candidate.is_absolute() else ROOT / candidate


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(BRIEF.is_file() and SOURCE.is_file(), "pending brief or source manifest missing")
    brief_raw = BRIEF.read_bytes()
    brief = _load_json_bytes(brief_raw, "pending brief")
    source_raw = SOURCE.read_bytes()
    source = _load_json_bytes(source_raw, "source manifest")
    _require(brief.get("component") == NAME, "brief component guard")
    _require(brief.get("phase") == PHASE, "brief phase guard")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief goal/one-entry guard",
    )
    checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    _require(_resolve_external_path(checkpoint.get("path")).resolve() == SOURCE.resolve(), "checkpoint path")
    _require(checkpoint.get("sha256") == _sha256_bytes(source_raw), "checkpoint hash")
    preimage = _resolve_external_path(checkpoint.get("preimage_path"))
    _require(preimage.is_file(), "brief preimage missing")
    _require(
        checkpoint.get("preimage_sha256") == _sha256_file(preimage),
        "brief preimage hash",
    )
    _require(source.get("component") == NAME, "manifest component guard")
    _require(source.get("one_entry_only") is True, "manifest one-entry guard")
    _require(source.get("previous_component") == PREVIOUS, "manifest predecessor guard")
    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "manifest input contract missing")
    _require(contract.get("leaf_api") == "build(c, *, response, cells, mu, joint, pulse)", "leaf API contract")
    _require(contract.get("fixture_api") == "run_fixture(leaf, c)", "fixture API contract")
    _require(contract.get("interval_precision_dps") == INTERVAL_DPS, "interval precision contract")
    policy = source.get("source_policy")
    _require(isinstance(policy, Mapping), "manifest source policy missing")
    for flag in (
        "authenticate_all_input_bytes_before_import_or_decode",
        "authenticate_raw_gzip_before_decode",
        "fixture_before_saved_receipt_decode",
        "exclusive_entry_and_result",
        "complete_result_and_failures_are_saved",
    ):
        _require(policy.get(flag) is True, "source policy guard: " + flag)
    return source, _manifest_hashes(source), _sha256_bytes(source_raw), _sha256_bytes(brief_raw)


def _all_authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    """Read and authenticate every input before any archive is decompressed."""

    raw: dict[str, bytes] = {}
    for key in FILES:
        data = (ROOT / key).read_bytes()
        _require(_sha256_bytes(data).lower() == hashes[key].lower(), "input changed after entry: " + key)
        raw[key] = data
    return raw


def _canonical(value: Any) -> bytes:
    return _json_bytes(_json_safe(value))


def _hash_pin(table: Any, path: str, expected: str, label: str) -> None:
    _require(isinstance(table, Mapping), label + " input_hashes missing")
    found = []
    for key in (path, Path(path).name):
        if key in table:
            found.append(table[key])
    _require(found, label + " does not pin " + path)
    for value in found:
        _require(
            isinstance(value, str)
            and SHA256_RE.fullmatch(value) is not None
            and value.lower() == expected.lower(),
            label + " has a changed pin for " + path,
        )


def _verify_source_pair(
    report: Mapping[str, Any],
    receipt: Mapping[str, Any],
    archive: str,
    check: str,
    hashes: Mapping[str, str],
    gate: str,
) -> Any:
    family = report.get("source_family")
    _require(family is not None, archive + " source_family missing")
    _require(
        receipt.get("source_family") is not None
        and _canonical(receipt.get("source_family")) == _canonical(family),
        check + " source_family mismatch",
    )
    _require(receipt.get("all_passed") is True, check + " all_passed is not true")
    _require(receipt.get(gate) is True, check + " gate is not true")
    if gate in report:
        _require(report[gate] is False, archive + " raw source scope was overclaimed")
    _hash_pin(receipt.get("input_hashes"), archive, hashes[archive], check)
    checks = report.get("checks")
    if checks is not None:
        _require(isinstance(checks, Mapping) and checks.get("passed") is True, archive + " checks")
    return family


def _verify_partial(report: Mapping[str, Any], hashes: Mapping[str, str]) -> tuple[Any, Mapping[str, Any]]:
    family = report.get("source_family")
    _require(family is not None, REPAIRED + " source_family missing")
    checks = report.get("checks")
    _require(isinstance(checks, Mapping) and checks.get("passed") is True, REPAIRED + " checks")
    pins = report.get("input_hashes")
    for path in (COEFF, COEFF_CHECK, FULL, FULL_CHECK, LINEAR, LINEAR_CHECK):
        _hash_pin(pins, path, hashes[path], REPAIRED)
    rows = report.get("partial_rows")
    _require(isinstance(rows, list), REPAIRED + " partial_rows missing")
    selected = [row for row in rows if isinstance(row, Mapping) and row.get("x") == "1.27"]
    _require(len(selected) == 1, REPAIRED + " requires one partial row x=1.27")
    row = selected[0]
    _require(
        row.get("source_family") is not None
        and _canonical(row.get("source_family")) == _canonical(family),
        REPAIRED + " selected row source_family mismatch",
    )
    _require(
        isinstance(row.get("angular_inverse_response"), list)
        and len(row["angular_inverse_response"]) == 3,
        REPAIRED + " angular response shape",
    )
    return family, row


def _verify_joint(report: Mapping[str, Any], hashes: Mapping[str, str]) -> tuple[Any, Mapping[str, Any]]:
    family = report.get("source_family")
    _require(family is not None, JOINT + " source_family missing")
    checks = report.get("checks")
    _require(isinstance(checks, Mapping) and checks.get("passed") is True, JOINT + " checks")
    cached = report.get("cached_component_source")
    _require(isinstance(cached, Mapping), JOINT + " cached source missing")
    _require(cached.get("report_sha256") == hashes[FULL], JOINT + " uses a different FULL hash")
    block = report.get("joint_reserved_minus_target_source")
    _require(isinstance(block, Mapping), JOINT + " joint source missing")
    _require(block.get("source_identity_before_boxing") is True, JOINT + " source identity gate")
    _require(block.get("independent_target_consumed_once") is True, JOINT + " target gate")
    _require(block.get("ordinary_Z_orders") == [0, 1], JOINT + " ordinary Z orders")
    _require(block.get("exact_response") == [0, 0, 1], JOINT + " exact response")
    normalization = block.get("common_normalization")
    _require(isinstance(normalization, Mapping) and normalization.get("U_Z_exact_zero") is True, JOINT + " normalization")
    coefficients = block.get("coefficient_C0_C1")
    _require(isinstance(coefficients, list) and len(coefficients) == 2, JOINT + " coefficient shape")
    return family, block


def _tuple4(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 4, label + " must be an exact mpf tuple")
    _require(all(type(item) is int for item in value), label + " has non-integer tuple fields")
    sign, mantissa, exponent, bits = tuple(value)
    _require(sign in (0, 1) and mantissa >= 0 and type(exponent) is int and bits >= 0, label + " tuple fields")
    if mantissa == 0:
        _require(sign == 0 and bits == 0, label + " noncanonical zero")
    else:
        _require(mantissa.bit_length() == bits, label + " bitcount")
    return sign, mantissa, exponent, bits


def _raw_interval_tuples(value: Any, label: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if isinstance(value, Mapping):
        if "source_enclosure_exact_mpf_tuples" in value:
            value = value["source_enclosure_exact_mpf_tuples"]
        elif "exact_mpf_tuples" in value:
            value = value["exact_mpf_tuples"]
        elif "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            value = (value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"])
        else:
            raise RunnerInputError(label + " has no exact endpoint tuples")
    _require(isinstance(value, (list, tuple)) and len(value) == 2, label + " interval shape")
    return _tuple4(value[0], label + "[0]"), _tuple4(value[1], label + "[1]")


def _singleton(c: Any, endpoint: tuple[int, ...]) -> Any:
    return c.make_mpf((tuple(endpoint), tuple(endpoint)))


def _make_interval(c: Any, value: Any, label: str) -> Any:
    lower, upper = _raw_interval_tuples(value, label)
    _require(_singleton(c, lower) <= _singleton(c, upper), label + " endpoints reversed")
    return c.make_mpf((lower, upper))


def _endpoint_values(c: Any, value: Any) -> tuple[Any, Any]:
    _require(hasattr(value, "_mpi_"), "comparison operand is not an interval")
    lower, upper = value._mpi_
    return _singleton(c, tuple(lower)), _singleton(c, tuple(upper))


def _interval_packet(c: Any, value: Any, origin: str) -> dict[str, Any]:
    _require(hasattr(value, "_mpi_"), origin + " is not an interval")
    lower, upper = value._mpi_
    return {
        "kind": "interval_enclosure",
        "origin": origin,
        "exact_mpf_tuples": [list(lower), list(upper)],
        "source_value_used": False,
        "midpoint_used_as_function_value": False,
        "old_error_used_as_function_value": False,
    }


def _interval_from_singletons(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((tuple(lower._mpi_[0]), tuple(upper._mpi_[0])))


def _lt(left: Any, right: Any) -> bool:
    return bool(left < right)


def _coverage_report(c: Any, records: list[Mapping[str, Any]]) -> dict[str, Any]:
    """Check coverage of open ``(-1, 1)`` without changing cell identity/order."""

    domain_lower = _singleton(c, (1, 1, 0, 1))
    domain_upper = _singleton(c, (0, 1, 0, 1))
    endpoints = []
    for record in records:
        lower, upper = _endpoint_values(c, record["clipped_Z"])
        endpoints.append({"record": record, "lower": lower, "upper": upper})

    ordered: list[int] = []
    for index, current in enumerate(endpoints):
        position = len(ordered)
        while position and _lt(current["lower"], endpoints[ordered[position - 1]]["lower"]):
            position -= 1
        ordered.insert(position, index)

    cursor = domain_lower
    gaps: list[Any] = []
    contributors: list[int] = []
    for index in ordered:
        entry = endpoints[index]
        if not _lt(cursor, entry["upper"]):
            continue
        if _lt(cursor, entry["lower"]):
            gap_upper = entry["lower"] if _lt(entry["lower"], domain_upper) else domain_upper
            if _lt(cursor, gap_upper):
                gaps.append(_interval_from_singletons(c, cursor, gap_upper))
            cursor = domain_upper
            break
        contributors.append(index)
        if _lt(cursor, entry["upper"]):
            cursor = entry["upper"]
        if not _lt(cursor, domain_upper):
            break
    if _lt(cursor, domain_upper):
        gaps.append(_interval_from_singletons(c, cursor, domain_upper))

    ordered_packets = []
    for index in ordered:
        record = endpoints[index]["record"]
        ordered_packets.append(
            {
                "source_cell_index": record["source_cell_index"],
                "cell_id": record["cell_id"],
                "clipped_Z": _interval_packet(c, record["clipped_Z"], "coverage cell clipped_Z"),
            }
        )
    return {
        "domain": "(-1,1)",
        "boundary_equal_is_allowed": True,
        "covered_open_domain": not gaps,
        "coverage_gap_count": len(gaps),
        "coverage_gaps": [
            _interval_packet(c, gap, "uncovered whole-Z gap") for gap in gaps
        ],
        "coverage_order_with_ids": ordered_packets,
        "contributing_source_cell_indices": [records[i]["source_cell_index"] for i in contributors],
        "cells_preserved_in_original_order": True,
        "function_points_materialized": False,
        "live_owner_claimed": False,
        "source_identity_proved": False,
    }


def _import_module(path: str, module_name: str) -> Any:
    _require(path in (LEAF, FIXTURE), "unexpected module import: " + path)
    spec = importlib.util.spec_from_file_location(module_name, ROOT / path)
    _require(spec is not None and spec.loader is not None, "module load failed: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _decode_archive(raw: bytes, label: str) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError(label + " is not a valid gzip archive") from exc
    return _load_json_bytes(decoded, label + " decoded JSON")


def _error(exc: BaseException) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def _consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    progress: dict[str, Any],
) -> dict[str, Any]:
    raw = _all_authenticated_bytes(hashes)
    progress["all_raw_inputs_authenticated_before_decode"] = True

    from mpmath.ctx_iv import MPIntervalContext

    c = MPIntervalContext()
    c.dps = INTERVAL_DPS
    leaf = _import_module(LEAF, NAME + "_leaf")
    fixture = _import_module(FIXTURE, NAME + "_fixture")
    fixture_runner = getattr(fixture, "run_fixture", None)
    _require(callable(fixture_runner), "fixture must export run_fixture")
    try:
        fixture_result = fixture_runner(leaf, c)
        progress["fixture"] = _json_safe(fixture_result, context=c, origin="fixture")
        _require(
            isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True,
            "synthetic fixture failed",
        )
    except BaseException as exc:
        progress["fixture_error"] = _error(exc)
        raise
    progress["fixture_completed_before_receipt_decode"] = True

    coverage_complete = _coverage_report(c, [
        {"source_cell_index": 0, "cell_id": "left", "clipped_Z": c.mpf([-1, 0])},
        {"source_cell_index": 1, "cell_id": "right", "clipped_Z": c.mpf([0, 1])},
    ])
    coverage_gap = _coverage_report(c, [
        {"source_cell_index": 0, "cell_id": "left", "clipped_Z": c.mpf([-1, "-0.5"])},
        {"source_cell_index": 1, "cell_id": "right", "clipped_Z": c.mpf([0, 1])},
    ])
    progress["coverage_fixture"] = {"complete": coverage_complete, "gap": coverage_gap}
    _require(coverage_complete["covered_open_domain"] is True and coverage_gap["covered_open_domain"] is False, "coverage fixture failed")

    progress["archive_decode_started"] = True
    partial = _load_json_bytes(raw[REPAIRED], REPAIRED)
    coefficient = _decode_archive(raw[COEFF], COEFF)
    coefficient_check = _load_json_bytes(raw[COEFF_CHECK], COEFF_CHECK)
    full_five = _decode_archive(raw[FULL], FULL)
    full_five_check = _load_json_bytes(raw[FULL_CHECK], FULL_CHECK)
    linear = _decode_archive(raw[LINEAR], LINEAR)
    linear_check = _load_json_bytes(raw[LINEAR_CHECK], LINEAR_CHECK)
    joint = _load_json_bytes(raw[JOINT], JOINT)

    partial_family, partial_row = _verify_partial(partial, hashes)
    coefficient_family = _verify_source_pair(
        coefficient,
        coefficient_check,
        COEFF,
        COEFF_CHECK,
        hashes,
        "current_original_n1_actual_reserved_coefficient_C1_numeric_enclosures_installed",
    )
    full_family = _verify_source_pair(
        full_five,
        full_five_check,
        FULL,
        FULL_CHECK,
        hashes,
        "current_original_true_complete_five_defect_C1_numeric_bridge_installed",
    )
    linear_family = _verify_source_pair(
        linear,
        linear_check,
        LINEAR,
        LINEAR_CHECK,
        hashes,
        "current_original_n1_reserved_normalized_linear_repair_operator_installed",
    )
    joint_family, joint_block = _verify_joint(joint, hashes)
    for label, family in (
        ("coefficient", coefficient_family),
        ("full_five", full_family),
        ("linear", linear_family),
        ("joint", joint_family),
    ):
        _require(_canonical(family) == _canonical(partial_family), label + " source_family differs")

    response = tuple(
        _make_interval(c, value["source_enclosure_exact_mpf_tuples"], f"partial angular response[{index}]")
        for index, value in enumerate(partial_row["angular_inverse_response"])
    )
    matrix = linear.get("actual_normalized_matrix_enclosures")
    _require(isinstance(matrix, Mapping) and "mu" in matrix, "linear mu missing")
    mu = _make_interval(c, matrix["mu"], "linear.actual_normalized_matrix_enclosures.mu")

    coefficient_sample = coefficient.get("samples", {}).get("whole_Z")
    _require(isinstance(coefficient_sample, Mapping), "coefficient whole_Z sample missing")
    cells_raw = coefficient_sample.get("cells")
    _require(isinstance(cells_raw, list) and cells_raw, "coefficient whole_Z cells missing")
    cells: list[dict[str, Any]] = []
    coverage_records: list[dict[str, Any]] = []
    for index, raw_cell in enumerate(cells_raw):
        _require(isinstance(raw_cell, Mapping), f"coefficient whole_Z cells[{index}] is not an object")
        _require("clipped_Z" in raw_cell and "normalized_theta_inputs" in raw_cell, f"cell[{index}] missing fields")
        clipped = _make_interval(c, raw_cell["clipped_Z"], f"coefficient whole_Z cells[{index}].clipped_Z")
        theta_raw = raw_cell["normalized_theta_inputs"]
        _require(isinstance(theta_raw, (list, tuple)) and len(theta_raw) >= 2, f"cell[{index}] theta orders")
        theta = tuple(
            tuple(
                _make_interval(c, theta_raw[order][component], f"cell[{index}] theta[{order}][{component}]")
                for component in range(2)
            )
            for order in range(2)
        )
        cell_id = raw_cell.get("cell_id", raw_cell.get("original_source_Z_cell"))
        if cell_id is None:
            cell_id = raw_cell.get("source_id", raw_cell.get("id"))
        _require(cell_id is not None, f"cell[{index}] has no stable id")
        payload = {
            "source_cell_index": index,
            "cell_id": cell_id,
            "source_id": cell_id,
            "source_path": [COEFF, "samples", "whole_Z", "cells", index],
            "source_packet": raw_cell,
            "uncertainty": True,
            "source_error_classes": ["saved_whole_Z_cell_enclosure"],
            "clipped_Z": clipped,
            "normalized_theta_inputs": theta,
            "raw_source_cell": raw_cell,
            "raw_source_record_preserved": True,
            "raw_source_path": raw_cell.get("source_path", raw_cell.get("path")),
            "raw_source_id": raw_cell.get("source_id", raw_cell.get("id", cell_id)),
            "raw_errors": raw_cell.get("errors"),
        }
        cells.append(payload)
        coverage_records.append(
            {"source_cell_index": index, "cell_id": cell_id, "clipped_Z": clipped}
        )
    coverage = _coverage_report(c, coverage_records)

    full_sample = full_five.get("samples", {}).get("whole_Z")
    _require(isinstance(full_sample, Mapping), "full-five whole_Z sample missing")
    pressure_components = full_sample.get("pressure_component_enclosures")
    _require(isinstance(pressure_components, Mapping), "full-five pressure components missing")
    m15 = pressure_components.get("m15")
    _require(isinstance(m15, Mapping), "full-five m15 pressure missing")
    pulse_raw = m15.get("pulse")
    _require(isinstance(pulse_raw, list) and len(pulse_raw) == 2, "full-five pulse must be a plain Z0/Z1 list")
    pulse = [
        _make_interval(c, value, f"full-five pressure pulse[{index}]")
        for index, value in enumerate(pulse_raw)
    ]
    joint_coefficients = joint_block["coefficient_C0_C1"]
    joint_values = tuple(
        _make_interval(c, value, f"joint J=-N coefficient_C0_C1[{index}]")
        for index, value in enumerate(joint_coefficients)
    )

    ledger = {
        "partial_row_source": REPAIRED,
        "partial_row_x": "1.27",
        "response_source_has_explicit_Z": False,
        "mu_independence_is_prior_contract": True,
        "coefficient_cells_source": COEFF,
        "coefficient_cell_count": len(cells),
        "coefficient_cell_indices_preserved": [cell["source_cell_index"] for cell in cells],
        "full_pressure_source": FULL,
        "joint_source": JOINT,
        "joint_full_source_hash": hashes[FULL],
        "P1_consumed": False,
        "u_Z_exact_zero": True,
        "function_points_materialized": False,
        "live_owner_claimed": False,
    }

    builder = getattr(leaf, "build", None)
    _require(callable(builder), "new leaf must export build")
    progress["actual_build_calls"] = 1
    built = builder(
        c,
        response={"response": response, "Z_independent": True, "x": "1.27", "source_path": [REPAIRED, "partial_rows", "x=1.27", "angular_inverse_response"]},
        cells=cells,
        mu=mu,
        joint=joint_values,
        pulse=pulse,
    )
    _require(isinstance(built, Mapping), "leaf.build result must be a mapping")
    leaf_claim_conflict = (
        not coverage["covered_open_domain"]
        and any(built.get(key) is True for key in ("whole_Z_claimed", "whole_Z_coverage_proved"))
    )

    return {
        "schema": NAME + "_result_v1",
        "precision_dps": INTERVAL_DPS,
        "source_family": partial_family,
        "input_hashes": dict(hashes),
        "selected_partial_row": _json_safe(partial_row, context=c, origin="saved partial row"),
        "response": [
            _interval_packet(c, value, "saved angular inverse response") for value in response
        ],
        "cells": _json_safe(cells, context=c, origin="authenticated whole-Z cell"),
        "mu": _interval_packet(c, mu, "linear actual_normalized_matrix_enclosures.mu"),
        "joint": [
            _interval_packet(c, value, "joint J=-N coefficient C0/C1") for value in joint_values
        ],
        "pulse": [
            _interval_packet(c, value, "full-five signed pressure pulse Z0/Z1") for value in pulse
        ],
        "coverage": coverage,
        "fixed_parameters": {
            "u_Z": 0,
            "ordinary_Z_orders": [0, 1],
            "P1_consumed": False,
            "response_query_x": "1.27",
            "response_Z_explicit": False,
            "mu_independence_is_prior_contract": True,
        },
        "source_lineage": {
            "partial": REPAIRED,
            "coefficient": COEFF,
            "coefficient_check": COEFF_CHECK,
            "full": FULL,
            "full_check": FULL_CHECK,
            "linear": LINEAR,
            "linear_check": LINEAR_CHECK,
            "joint": JOINT,
            "P_and_J_share_full_hash": hashes[FULL],
        },
        "error_ledger": ledger,
        "build": _json_safe(built, context=c, origin="leaf build result"),
        "actual_build_calls": 1,
        "whole_Z_input_cell_coverage_proved": bool(coverage["covered_open_domain"]),
        "whole_Z_claimed": False,
        "leaf_whole_Z_claim_conflict": leaf_claim_conflict,
        "function_points_materialized": False,
        "live_owner_claimed": False,
        "source_identity_proved": False,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": NAME + "_entry_v1",
            "component": NAME,
            "phase": PHASE,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "entered": True,
            "all_raw_inputs_authenticated_before_decode": True,
            "fixture_before_receipt_decode": True,
            "actual_build_calls": 0,
            "function_points_materialized": False,
            "live_owner_claimed": False,
        },
    )
    out: dict[str, Any] = {
        "schema": NAME + "_result_v1",
        "component": NAME,
        "phase": PHASE,
        "source_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "source_manifest": source,
        "status": "failed",
        "archive_decode_started": False,
        "all_raw_inputs_authenticated_before_decode": True,
        "fixture_before_receipt_decode": True,
        "actual_build_calls": 0,
        "fixture": None,
        "function_points_materialized": False,
        "live_owner_claimed": False,
    }
    try:
        output = _consume(source, hashes, out)
    except BaseException as exc:
        out["error"] = _error(exc)
    else:
        out["status"] = "returned"
        out["output"] = output
    result_sha = _write_exclusive(RESULT, out)
    saved_error = out.get("error")
    error_message = saved_error.get("message") if isinstance(saved_error, Mapping) else None
    print(
        json.dumps(
            {
                "component": NAME,
                "phase": PHASE,
                "status": out["status"],
                "result_sha256": result_sha,
                "error": error_message,
            },
            sort_keys=True,
        )
    )
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())


