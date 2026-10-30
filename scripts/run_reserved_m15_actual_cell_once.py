"""Freeze and consume the saved reserved-m15 source cells exactly once.

This runner is the persistence boundary for the saved source-cell assembly.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "reserved_m15_actual_cell_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF = Path("src/openai_ns_reconstruction/reserved_m15_cell_assembly.py")
RUNNER = Path("scripts/run_reserved_m15_actual_cell_once.py")
PARTIAL = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_repaired_Rp_m15_numeric.json"
)
COEFFICIENT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_coefficient_enclosures.json.gz"
)
COEFFICIENT_CHECK = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_coefficient_enclosures_check.json"
)
FULL_FIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge.json.gz"
)
FULL_FIVE_CHECK = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge_check.json"
)
LINEAR = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator.json.gz"
)
LINEAR_CHECK = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator_check.json"
)
REQUIRED_INPUTS = frozenset(
    {
        LEAF.as_posix(), RUNNER.as_posix(), PARTIAL.as_posix(),
        COEFFICIENT.as_posix(), COEFFICIENT_CHECK.as_posix(),
        FULL_FIVE.as_posix(), FULL_FIVE_CHECK.as_posix(),
        LINEAR.as_posix(), LINEAR_CHECK.as_posix(),
    }
)
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
INTERVAL_DPS = 130
FIXED_SHA256 = {
    PARTIAL.as_posix(): "bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871",
}


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
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


def _load_json_file(path: Path, label: str) -> dict[str, Any]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise RunnerInputError(label + " cannot be read") from exc
    return _load_json_bytes(raw, label)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value, ensure_ascii=False, sort_keys=True,
            separators=(",", ":"), allow_nan=False,
        ).encode("utf-8") + b"\n"
    )


def _json_safe(value: Any, *, context: Any = None) -> Any:
    if context is not None and hasattr(value, "_mpi_"):
        return _interval_packet(context, value, "directed interval")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item, context=context) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item, context=context) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item, context=context) for item in value), key=repr)
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
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(not relative.is_absolute() and ".." not in relative.parts,
             label + " must be repository-relative")
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(),
             label + " is missing or escapes repository")
    return key, path


def _manifest_records(source: Mapping[str, Any]) -> list[tuple[Any, Any, str]]:
    records: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            records.append((item.get("path"), item.get("sha256"), f"source.files[{index}]"))
    for field in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            records.extend((key, value, f"source.{field}") for key, value in table.items())
    _require(records, "source manifest contains no raw SHA entries")
    return records


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for raw_path, expected, label in _manifest_records(source):
        key, path = _repository_path(raw_path, label)
        _require(isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
                 "invalid SHA-256 for " + key)
        digest = expected.lower()
        fixed = FIXED_SHA256.get(key)
        if fixed is not None:
            _require(digest == fixed, "fixed source SHA pin changed: " + key)
        if key in hashes:
            _require(hashes[key].lower() == digest, "conflicting SHA entries for " + key)
            continue
        actual = _sha256_file(path)
        _require(actual.lower() == digest, "raw SHA changed: " + key)
        hashes[key] = actual
    missing = sorted(REQUIRED_INPUTS - set(hashes))
    _require(not missing, "manifest omits frozen input(s): " + ", ".join(missing))
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
    _require(candidate.resolve() == SOURCE.resolve(), "source checkpoint points elsewhere")
    expected = checkpoint.get("sha256")
    actual = _sha256_file(SOURCE)
    _require(isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None
             and actual.lower() == expected.lower(), "source checkpoint SHA changed")
    return actual


def _guard_true(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "source manifest guard missing: " + label)
    _require(any(source[name] is True for name in present),
             "source manifest guard is not true: " + label)


def _guard_false_if_present(
    source: Mapping[str, Any], names: tuple[str, ...], label: str
) -> None:
    for name in names:
        if name in source:
            _require(source[name] is False, "source manifest guard changed: " + label)


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief = _load_json_file(BRIEF, "pending brief")
    source = _load_json_file(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "pending brief component changed")
    _require(brief.get("phase") == PHASE, "pending brief phase changed")
    _require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False
             and brief.get("do_not_reenter") is True, "pending brief goal flags changed")
    source_sha = _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "source manifest component changed")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")
    _guard_true(
        source,
        ("all_raw_inputs_authenticated_before_decompression",
         "all_raw_inputs_authenticated_before_decode",
         "raw_sha256_authenticated_before_entry",
         "all_input_raw_sha256_authenticated_before_decode"),
        "raw archive authentication",
    )
    _guard_true(
        source,
        ("no_old_producer_import", "no_old_helper_import", "no_ancestor_replay",
         "no_old_producer_or_partial_row_call", "no_old_partial_row_call",
         "no_old_producer_or_helper", "no_weighted_beta_prefix_import",
         "no_weighted_beta_prefix_call"),
        "legacy producer/helper exclusion",
    )
    _guard_true(
        source,
        ("protected_graph_nodes_not_read", "no_protected_graph_node_read",
         "protected_nodes_not_indexed", "no_protected_or_target_graph_read",
         "no_protected_graph_node_reads", "protected_graph_nodes_unread"),
        "protected graph exclusion",
    )
    _guard_false_if_present(
        source,
        ("graph_traversal", "graph_edge_search", "source_graph_traversed",
         "protected_graph_node_read", "protected_nodes_indexed"),
        "graph/protected-node scope",
    )
    hashes = _manifest_hashes(source)
    _require(LEAF.as_posix() in hashes, "new assembly leaf is not pinned")
    _require(RUNNER.as_posix() in hashes, "new runner is not pinned")
    return source, hashes, source_sha


def _authenticated_bytes(relative: Path, hashes: Mapping[str, str]) -> bytes:
    key = relative.as_posix()
    _require(key in hashes, "raw input is not pinned: " + key)
    raw = (ROOT / relative).read_bytes()
    _require(_sha256_bytes(raw).lower() == hashes[key].lower(),
             "raw input changed after entry: " + key)
    return raw


def _all_authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    """Read/hash every declared input before any gzip member is decoded."""

    return {
        key: _authenticated_bytes(Path(key), hashes)
        for key in sorted(hashes)
    }


def _canonical(value: Any) -> bytes:
    return _json_bytes(_json_safe(value))


def _hash_pin(table: Any, relative: Path, expected: str, label: str) -> None:
    _require(isinstance(table, Mapping), label + " input_hashes missing")
    keys = (relative.as_posix(), relative.name)
    found = [(key, table[key]) for key in keys if key in table]
    _require(found, label + " does not pin " + relative.as_posix())
    for key, value in found:
        _require(isinstance(value, str) and SHA256_RE.fullmatch(value) is not None
                 and value.lower() == expected.lower(),
                 label + " has a changed pin for " + key)


def _verify_source_pair(
    report: Mapping[str, Any],
    receipt: Mapping[str, Any],
    archive: Path,
    check: Path,
    hashes: Mapping[str, str],
    gate: str,
) -> Any:
    family = report.get("source_family")
    _require(family is not None, archive.as_posix() + " source_family missing")
    _require(receipt.get("source_family") is not None
             and _canonical(receipt.get("source_family")) == _canonical(family),
             check.as_posix() + " source_family mismatch")
    _require(receipt.get("all_passed") is True,
             check.as_posix() + " all_passed is not true")
    _require(receipt.get(gate) is True, check.as_posix() + " gate is not true")
    if gate in report:
        _require(report[gate] is False,
                 archive.as_posix() + " raw source scope was overclaimed")
    _hash_pin(receipt.get("input_hashes"), archive, hashes[archive.as_posix()],
              check.as_posix())
    checks = report.get("checks")
    if checks is not None:
        _require(isinstance(checks, Mapping) and checks.get("passed") is True,
                 archive.as_posix() + " checks.passed is not true")
    return family


def _verify_partial(report: Mapping[str, Any], hashes: Mapping[str, str]) -> Any:
    family = report.get("source_family")
    _require(family is not None, PARTIAL.as_posix() + " source_family missing")
    checks = report.get("checks")
    _require(isinstance(checks, Mapping) and checks.get("passed") is True,
             PARTIAL.as_posix() + " checks.passed is not true")
    pins = report.get("input_hashes")
    for path in (COEFFICIENT, COEFFICIENT_CHECK, FULL_FIVE, FULL_FIVE_CHECK,
                 LINEAR, LINEAR_CHECK):
        _hash_pin(pins, path, hashes[path.as_posix()], PARTIAL.as_posix())
    rows = report.get("partial_rows")
    _require(isinstance(rows, list), PARTIAL.as_posix() + " partial_rows missing")
    selected = [row for row in rows
                if isinstance(row, Mapping) and row.get("x") == "1.27"]
    _require(len(selected) == 1,
             PARTIAL.as_posix() + " requires one partial row with x='1.27'")
    row = selected[0]
    _require(row.get("source_family") is not None
             and _canonical(row.get("source_family")) == _canonical(family),
             PARTIAL.as_posix() + " selected row source_family mismatch")
    _require(isinstance(row.get("angular_inverse_response"), list)
             and len(row["angular_inverse_response"]) == 3,
             PARTIAL.as_posix() + " angular response shape changed")
    _require(isinstance(row.get("coefficient_C0_C1"), list)
             and len(row["coefficient_C0_C1"]) == 2,
             PARTIAL.as_posix() + " coefficient shape changed")
    return family, row


def _tuple4(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 4,
             label + " must be one exact mpf tuple")
    _require(all(type(item) is int for item in value),
             label + " has non-integer mpf tuple fields")
    sign, mantissa, exponent, bits = tuple(value)
    _require(sign in (0, 1) and mantissa >= 0 and type(exponent) is int and bits >= 0,
             label + " has invalid mpf tuple fields")
    if mantissa == 0:
        _require(sign == 0 and bits == 0,
                 label + " has a noncanonical zero mpf tuple")
    else:
        _require(mantissa.bit_length() == bits,
                 label + " has an invalid mpf bitcount")
    return sign, mantissa, exponent, bits


def _raw_interval_tuples(
    value: Any, label: str
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if isinstance(value, Mapping):
        if "source_enclosure_exact_mpf_tuples" in value:
            value = value["source_enclosure_exact_mpf_tuples"]
        elif "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            value = (value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"])
        else:
            raise RunnerInputError(label + " has no exact mpf interval tuples")
    _require(isinstance(value, (list, tuple)) and len(value) == 2,
             label + " must contain exactly two exact mpf endpoints")
    return _tuple4(value[0], label + "[0]"), _tuple4(value[1], label + "[1]")


def _make_interval(ctx: Any, mp: Any, value: Any, label: str) -> Any:
    lower, upper = _raw_interval_tuples(value, label)
    lower_value = mp.make_mpf(lower)
    upper_value = mp.make_mpf(upper)
    _require(lower_value <= upper_value, label + " endpoints are reversed")
    return ctx.make_mpf((lower, upper))


def _interval_packet(ctx: Any, value: Any, origin: str) -> dict[str, Any]:
    _require(hasattr(value, "_mpi_"), origin + " is not an interval value")
    lower, upper = value._mpi_
    import mpmath as mp
    return {
        "kind": "interval_enclosure",
        "origin": origin,
        "scientific": [
            mp.nstr(mp.make_mpf(tuple(lower)), 30, min_fixed=0, max_fixed=0),
            mp.nstr(mp.make_mpf(tuple(upper)), 30, min_fixed=0, max_fixed=0),
        ],
        "exact_mpf_tuples": [list(lower), list(upper)],
        "source_value_used": False,
        "midpoint_used_as_function_value": False,
        "old_error_used_as_function_value": False,
    }


def _endpoint_values(mp: Any, value: Any) -> tuple[Any, Any]:
    _require(hasattr(value, "_mpi_"), "comparison operand is not an interval")
    lower, upper = value._mpi_
    return mp.make_mpf(tuple(lower)), mp.make_mpf(tuple(upper))


def _directed_hull(ctx: Any, mp: Any, values: Sequence[Any]) -> Any:
    _require(values, "directed hull requires at least one interval")
    lowers = [value._mpi_[0] for value in values]
    uppers = [value._mpi_[1] for value in values]
    lower = min(lowers, key=lambda raw: mp.make_mpf(tuple(raw)))
    upper = max(uppers, key=lambda raw: mp.make_mpf(tuple(raw)))
    return ctx.make_mpf((tuple(lower), tuple(upper)))


def _width(ctx: Any, value: Any) -> Any:
    lower, upper = value._mpi_
    hi = ctx.make_mpf((tuple(upper), tuple(upper)))
    lo = ctx.make_mpf((tuple(lower), tuple(lower)))
    return hi - lo


def _import_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new assembly leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "reserved_m15_actual_cell_v1_new_leaf", path
    )
    _require(spec is not None and spec.loader is not None,
             "new assembly leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new assembly callable is missing: " + name)
    return function


def _decode_archive(raw: bytes, label: str) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError(label + " is not a valid gzip archive") from exc
    return _load_json_bytes(decoded, label + " decoded JSON")


def _consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    stats: dict[str, Any],
) -> dict[str, Any]:
    # Authenticate every declared raw input before any gzip member is decoded.
    # In particular, no exact_graph_nodes field is inspected.
    raw = _all_authenticated_bytes(hashes)
    stats["archive_decode_started"] = True
    partial = _load_json_bytes(raw[PARTIAL.as_posix()], PARTIAL.as_posix())
    coefficient = _decode_archive(
        raw[COEFFICIENT.as_posix()], COEFFICIENT.as_posix()
    )
    coefficient_check = _load_json_bytes(
        raw[COEFFICIENT_CHECK.as_posix()], COEFFICIENT_CHECK.as_posix()
    )
    full_five = _decode_archive(raw[FULL_FIVE.as_posix()], FULL_FIVE.as_posix())
    full_five_check = _load_json_bytes(
        raw[FULL_FIVE_CHECK.as_posix()], FULL_FIVE_CHECK.as_posix()
    )
    linear = _decode_archive(raw[LINEAR.as_posix()], LINEAR.as_posix())
    linear_check = _load_json_bytes(
        raw[LINEAR_CHECK.as_posix()], LINEAR_CHECK.as_posix()
    )

    partial_family, partial_row = _verify_partial(partial, hashes)
    coefficient_family = _verify_source_pair(
        coefficient, coefficient_check, COEFFICIENT, COEFFICIENT_CHECK, hashes,
        "current_original_n1_actual_reserved_coefficient_C1_numeric_enclosures_installed",
    )
    full_family = _verify_source_pair(
        full_five, full_five_check, FULL_FIVE, FULL_FIVE_CHECK, hashes,
        "current_original_true_complete_five_defect_C1_numeric_bridge_installed",
    )
    linear_family = _verify_source_pair(
        linear, linear_check, LINEAR, LINEAR_CHECK, hashes,
        "current_original_n1_reserved_normalized_linear_repair_operator_installed",
    )
    for name, family in (("coefficient", coefficient_family),
                         ("full_five", full_family), ("linear", linear_family)):
        _require(_canonical(family) == _canonical(partial_family),
                 name + " source_family differs from repaired-Rp report")

    coefficient_sample = coefficient.get("samples", {}).get("whole_Z")
    _require(isinstance(coefficient_sample, Mapping), "coefficient whole_Z sample missing")
    cells_raw = coefficient_sample.get("cells")
    _require(isinstance(cells_raw, list) and cells_raw,
             "coefficient whole_Z cells missing")
    full_sample = full_five.get("samples", {}).get("whole_Z")
    _require(isinstance(full_sample, Mapping), "full-five whole_Z sample missing")
    pressure_raw = full_sample.get(
        "normalized_full_pressure_moment_C1_enclosures", {}
    ).get("m15")
    _require(pressure_raw is not None, "full-five whole_Z m15 pressure pair missing")
    matrix_raw = linear.get("actual_normalized_matrix_enclosures")
    _require(isinstance(matrix_raw, Mapping) and "mu" in matrix_raw,
             "linear normalized matrix mu missing")

    from mpmath.ctx_iv import MPIntervalContext
    import mpmath as mp

    ctx = MPIntervalContext()
    ctx.dps = INTERVAL_DPS
    mu = _make_interval(ctx, mp, matrix_raw["mu"],
                        "linear.actual_normalized_matrix_enclosures.mu")
    amplitude = ctx.exp(-ctx.mpf(17) / 10 - ctx.mpf(5) * mu / 2)
    amplitude_z1 = ctx.mpf(0)

    response_values = []
    for index, row in enumerate(partial_row["angular_inverse_response"]):
        _require(isinstance(row, Mapping)
                 and "source_enclosure_exact_mpf_tuples" in row,
                 f"partial_rows[x=1.27].angular_inverse_response[{index}]"
                 " lacks exact tuples")
        response_values.append(_make_interval(
            ctx, mp, row["source_enclosure_exact_mpf_tuples"],
            f"partial_rows[x=1.27].angular_inverse_response[{index}]",
        ))
    response = tuple(response_values)
    pressure = tuple(
        _make_interval(ctx, mp, value, f"full-five whole_Z m15 pressure[{index}]")
        for index, value in enumerate(pressure_raw)
    )
    _require(len(pressure) == 2, "full-five m15 pressure pair is malformed")
    old_coefficient = tuple(
        _make_interval(ctx, mp, value,
                       f"partial_rows[x=1.27].coefficient_C0_C1[{index}]")
        for index, value in enumerate(partial_row["coefficient_C0_C1"])
    )

    source_cells: list[dict[str, Any]] = []
    selected_raw_cells: list[Any] = []
    cell_packets: list[dict[str, Any]] = []
    for index, raw_cell in enumerate(cells_raw):
        _require(isinstance(raw_cell, Mapping),
                 f"coefficient whole_Z cells[{index}] is not an object")
        _require("clipped_Z" in raw_cell and "normalized_theta_inputs" in raw_cell,
                 f"coefficient whole_Z cells[{index}] is incomplete")
        clipped = _make_interval(
            ctx, mp, raw_cell["clipped_Z"],
            f"coefficient whole_Z cells[{index}].clipped_Z",
        )
        clipped_lo, clipped_hi = _endpoint_values(mp, clipped)
        if not (clipped_lo <= mp.mpf(0) <= clipped_hi):
            continue
        theta_raw = raw_cell["normalized_theta_inputs"]
        _require(isinstance(theta_raw, (list, tuple)) and len(theta_raw) >= 2,
                 f"coefficient whole_Z cells[{index}].normalized_theta_inputs"
                 " is incomplete")
        theta = tuple(
            tuple(
                _make_interval(
                    ctx, mp, theta_raw[channel][component],
                    f"coefficient whole_Z cells[{index}].normalized_theta_inputs"
                    f"[{channel}][{component}]",
                )
                for component in range(2)
            )
            for channel in range(2)
        )
        cell_id = raw_cell.get("cell_id", raw_cell.get("original_source_Z_cell"))
        _require(cell_id is not None,
                 f"coefficient whole_Z cells[{index}] has no stable cell id")
        ledger = {"source_cell_index": index, "raw_source_cell": raw_cell}
        source_cells.append({
            "cell_id": cell_id,
            "clipped_Z": raw_cell["clipped_Z"],
            "normalized_theta_inputs": theta,
            "error_ledger": ledger,
        })
        selected_raw_cells.append(raw_cell)
        cell_packets.append({
            "cell_id": _json_safe(cell_id),
            "clipped_Z": _interval_packet(
                ctx, clipped, "authenticated coefficient cell clipped_Z"
            ),
            "normalized_theta_inputs": [
                [
                    _interval_packet(
                        ctx, theta[channel][component],
                        "authenticated coefficient normalized theta interval",
                    )
                    for component in range(2)
                ]
                for channel in range(2)
            ],
            "ledger_preserved": True,
        })
    _require(source_cells, "no coefficient whole_Z cell covers Z=0")

    ledger = {
        "partial_row": partial_row,
        "coefficient_whole_Z_selected_cells": selected_raw_cells,
        "full_five_whole_Z_pressure_source": pressure_raw,
        "linear_mu_source": matrix_raw["mu"],
        "old_212_term_state": "untouched; no legacy ledger was interpreted or rewritten",
    }
    assembler = _callable(_import_leaf(), "assemble_source_cells")
    hull_calls = 0

    def hull(values: Sequence[Any]) -> Any:
        nonlocal hull_calls
        hull_calls += 1
        return _directed_hull(ctx, mp, values)

    stats["assembly_invocations"] = 1
    assembled = assembler(
        response=response,
        source_cells=source_cells,
        a=amplitude,
        ordinary_z1=amplitude_z1,
        m15_pressure=pressure,
        hull=hull,
        error_ledger=ledger,
    )
    _require(hull_calls == 2,
             "assembly did not request one directed hull per C0/C1 component")
    new_hull = tuple(assembled.get("hull_C0_C1", ()))
    _require(len(new_hull) == 2 and all(hasattr(value, "_mpi_") for value in new_hull),
             "assembly returned malformed hull_C0_C1")

    comparison = []
    for index, (old_value, new_value) in enumerate(zip(old_coefficient, new_hull)):
        old_lo, old_hi = _endpoint_values(mp, old_value)
        new_lo, new_hi = _endpoint_values(mp, new_value)
        contained = bool(new_lo >= old_lo and new_hi <= old_hi)
        old_width = _width(ctx, old_value)
        new_width = _width(ctx, new_value)
        comparison.append({
            "component": "C0" if index == 0 else "C1",
            "old_saved_interval": _interval_packet(
                ctx, old_value, "saved partial-row coefficient interval"
            ),
            "new_directed_hull": _interval_packet(
                ctx, new_value, "new source-cell directed hull"
            ),
            "new_enclosure_within_old": contained,
            "old_width": _interval_packet(ctx, old_width, "saved interval width"),
            "new_width": _interval_packet(ctx, new_width, "new hull interval width"),
            "width_delta_new_minus_old": _interval_packet(
                ctx, new_width - old_width, "interval width difference"
            ),
            "possible_cross_precision_rounding_noncontainment": not contained,
            "comparison_is_interval_width_evidence_only": True,
            "source_identity_proved": False,
        })

    return {
        "precision_dps": INTERVAL_DPS,
        "selected_partial_row": partial_row,
        "legacy_partial_row_raw_untouched": partial_row,
        "legacy_partial_row_approximations_not_used": True,
        "selected_partial_row_scope": {
            "x": "1.27",
            "query_Z": "0",
            "scope_basis": "check_partial_rows source contract; saved row was not recomputed",
            "whole_Z_sample_not_used_for_partial_row": True,
        },
        "source_family": partial_family,
        "authenticated_input_intervals": {
            "angular_inverse_response": [
                _interval_packet(ctx, value, "authenticated saved angular response interval")
                for value in response
            ],
            "normalized_theta_inputs_by_cell": cell_packets,
            "m15_pressure": [
                _interval_packet(ctx, value, "authenticated full-five m15 pressure interval")
                for value in pressure
            ],
            "mu": _interval_packet(ctx, mu, "authenticated linear-operator mu interval"),
        },
        "computed_directed_intervals": {
            "A": _interval_packet(
                ctx, amplitude,
                "A=exp(-17/10-5*mu/2) from authenticated mu interval",
            ),
            "A_Z": _interval_packet(ctx, amplitude_z1, "exact Z=0 amplitude derivative"),
        },
        "source_cells": cell_packets,
        "error_ledger": ledger,
        "assembly": _json_safe(assembled, context=ctx),
        "comparison": comparison,
        "hull_calls": hull_calls,
        "legacy_ledger_processed": False,
        "legacy_212_term_state_untouched": True,
        "midpoint_used_as_function_value": False,
        "old_error_used_as_function_value": False,
        "source_identity_proved": False,
        "ordinary_z_source_identity_proved": False,
        "protected_graph_node_read": False,
        "old_producer_or_partial_helper_called": False,
    }


def main() -> int:
    source, hashes, source_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "manifest_input_count": len(hashes),
            "archive_decode_started": False,
            "assembly_invocations": 0,
            "midpoint_used_as_function_value": False,
            "old_error_used_as_function_value": False,
            "protected_graph_node_read": False,
        },
    )
    stats = {"archive_decode_started": False, "assembly_invocations": 0}
    status = "failed"
    packet: Any = None
    error: Any = None
    try:
        packet = _consume(source, hashes, stats)
        status = "returned"
    except BaseException as exc:
        error = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    output = {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": status,
        "source_checkpoint_sha256": source_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "manifest_input_count": len(hashes),
        "result": packet,
        "error": error,
        "source_manifest": source,
        "archive_decode_started": stats["archive_decode_started"],
        "assembly_invocations": stats["assembly_invocations"],
        "midpoint_used_as_function_value": False,
        "old_error_used_as_function_value": False,
        "protected_graph_node_read": False,
        "old_producer_or_partial_helper_called": False,
    }
    result_sha = _write_exclusive(RESULT, output)
    print(json.dumps(
        {
            "component": COMPONENT,
            "phase": PHASE,
            "status": status,
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "result_sha256": result_sha,
            "input_hashes": hashes,
        },
        sort_keys=True,
    ))
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
