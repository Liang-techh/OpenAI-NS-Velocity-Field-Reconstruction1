"""Run the variable-Z reserved response consumer exactly once.

The prior v3 result is an authenticated input packet, not a function value.
The independent fixture runs before that JSON is decoded.  Each saved source
cell is passed to the new leaf exactly once; unknown ``r_i_Z`` and ``mu_Z``
terms are retained as explicit sensitivity coefficients and are not set to
zero as a mathematical assumption.
"""

from __future__ import annotations

import dataclasses
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
NAME = "reserved_m15_variable_z"
PHASE = NAME + "_pending"
PREVIOUS = "reserved_partial_wholez_v3"
DIRECTORY = ROOT / "tmp" / NAME
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF = "src/openai_ns_reconstruction/reserved_m15_variable_z.py"
FIXTURE = "src/openai_ns_reconstruction/reserved_m15_variable_z_fixture.py"
RUNNER = "scripts/run_reserved_variable_z_once.py"
PREVIOUS_RESULT = "tmp/reserved_partial_wholez_v3/result.json"
PERMANENT_RECEIPT = "experiments/root_st073/reserved_partial_wholez_v3_receipt.json.gz"
FILES = (LEAF, FIXTURE, RUNNER, PREVIOUS_RESULT, PERMANENT_RECEIPT)

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
INTERVAL_DPS = 130
FIXED_SHA256 = {
    PREVIOUS_RESULT: "5bd627d0744d36804547c6ebbf02884a68b4163ee4eb795ee29b35ecd561a2cb",
    PERMANENT_RECEIPT: "967fd2f9fe3fb8361b2bf632bac2987213c350056b1e3aa7ad5c83056cbf1308",
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


def _as_fields(value: Any, label: str) -> dict[str, Any]:
    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items()}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: getattr(value, field.name)
            for field in dataclasses.fields(value)
        }
    raise RunnerInputError(label + " is not a dataclass or mapping")


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
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _json_safe(
                getattr(value, field.name), context=context, origin=origin + "." + field.name
            )
            for field in dataclasses.fields(value)
        }
    if isinstance(value, Mapping):
        return {
            str(key): _json_safe(item, context=context, origin=origin + "." + str(key))
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
        raise RunnerInputError("one-entry receipt already exists: " + str(path)) from exc
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


def _resolve_external_path(raw_path: Any) -> Path:
    _require(isinstance(raw_path, str) and raw_path, "checkpoint path missing")
    candidate = Path(raw_path)
    return candidate if candidate.is_absolute() else ROOT / candidate


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
        if key in (PREVIOUS_RESULT, PERMANENT_RECEIPT):
            _require(
                row.get("raw_sha256_before_decode") is True,
                "raw input was not authenticated: " + key,
            )
        hashes[key] = actual
    _require(set(hashes) == set(FILES), "manifest input set mismatch")
    return hashes


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
    _require(
        _resolve_external_path(checkpoint.get("path")).resolve() == SOURCE.resolve(),
        "checkpoint path",
    )
    _require(checkpoint.get("sha256") == _sha256_bytes(source_raw), "checkpoint hash")
    preimage = _resolve_external_path(checkpoint.get("preimage_path"))
    _require(preimage.is_file(), "brief preimage missing")
    _require(checkpoint.get("preimage_sha256") == _sha256_file(preimage), "brief preimage hash")

    _require(source.get("component") == NAME, "manifest component guard")
    _require(source.get("previous_component") == PREVIOUS, "manifest predecessor guard")
    _require(source.get("one_entry_only") is True, "manifest one-entry guard")
    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "manifest input contract missing")
    _require(
        contract.get("leaf_api")
        == "evaluate_partial_reserved_minus_pulse(iv, r_i, r_i_Z, mu, mu_Z, theta_i, theta_i_Z, J, J_Z, P, P_Z, Z)",
        "leaf API contract",
    )
    _require(contract.get("fixture_api") == "run_fixture(leaf, c)", "fixture API contract")
    _require(contract.get("interval_precision_dps") == INTERVAL_DPS, "precision contract")
    _require(contract.get("previous_result_sha256") == FIXED_SHA256[PREVIOUS_RESULT], "previous result pin")
    _require(contract.get("permanent_receipt_raw_sha256") == FIXED_SHA256[PERMANENT_RECEIPT], "receipt pin")
    _require(contract.get("unknown_sensitivity_coefficients") == ["r0_Z", "r1_Z", "r2_Z", "mu_Z"], "sensitivity contract")
    policy = source.get("source_policy")
    _require(isinstance(policy, Mapping), "manifest source policy missing")
    for flag in (
        "authenticate_all_input_bytes_before_import_or_decode",
        "authenticate_previous_result_raw_before_decode",
        "authenticate_permanent_receipt_raw_before_decode",
        "fixture_before_input_json_decode",
        "exclusive_entry_and_result",
        "complete_result_and_failures_are_saved",
        "actual_leaf_call_once_per_cell",
    ):
        _require(policy.get(flag) is True, "source policy guard: " + flag)
    return source, _manifest_hashes(source), _sha256_bytes(source_raw), _sha256_bytes(brief_raw)


def _all_authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for key in FILES:
        data = (ROOT / key).read_bytes()
        _require(_sha256_bytes(data).lower() == hashes[key].lower(), "input changed after entry: " + key)
        raw[key] = data
    return raw


def _import_module(path: str, module_name: str) -> Any:
    _require(path in (LEAF, FIXTURE), "unexpected module import: " + path)
    spec = importlib.util.spec_from_file_location(module_name, ROOT / path)
    _require(spec is not None and spec.loader is not None, "module load failed: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _tuple4(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 4, label + " must be an exact mpf tuple")
    _require(all(type(item) is int for item in value), label + " has non-integer tuple fields")
    sign, mantissa, exponent, bits = tuple(value)
    _require(sign in (0, 1) and mantissa >= 0 and bits >= 0, label + " tuple fields")
    if mantissa == 0:
        _require(sign == 0 and bits == 0, label + " noncanonical zero")
    else:
        _require(mantissa.bit_length() == bits, label + " bitcount")
    return sign, mantissa, exponent, bits


def _raw_interval_tuples(value: Any, label: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if isinstance(value, Mapping):
        if "exact_mpf_tuples" in value:
            value = value["exact_mpf_tuples"]
        elif "source_enclosure_exact_mpf_tuples" in value:
            value = value["source_enclosure_exact_mpf_tuples"]
        elif "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            value = (value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"])
        else:
            raise RunnerInputError(label + " has no exact endpoint tuples")
    _require(isinstance(value, (list, tuple)) and len(value) == 2, label + " interval shape")
    return _tuple4(value[0], label + "[0]"), _tuple4(value[1], label + "[1]")


def _make_interval(c: Any, value: Any, label: str) -> Any:
    lower, upper = _raw_interval_tuples(value, label)
    interval = c.make_mpf((lower, upper))
    lower_value = c.make_mpf((lower, lower))
    upper_value = c.make_mpf((upper, upper))
    _require(bool(lower_value <= upper_value), label + " endpoints reversed")
    return interval


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


def _decode_archive(raw: bytes, label: str) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError(label + " is not a valid gzip archive") from exc
    return _load_json_bytes(decoded, label + " decoded JSON")


def _canonical(value: Any) -> bytes:
    return _json_bytes(value)


def _receipt_source_result(receipt: Mapping[str, Any], previous: Mapping[str, Any], previous_sha: str) -> None:
    """Authenticate the old archive shape and compare its saved result directly."""

    for key in ("component", "entry", "frozen_code_base64", "result", "source"):
        _require(key in receipt, "permanent receipt missing top-level " + key)
    _require(receipt.get("component") == PREVIOUS, "permanent receipt component mismatch")
    _require(receipt.get("result") is not None, "permanent receipt does not preserve result")
    candidate = receipt["result"]
    _require(isinstance(candidate, Mapping), "permanent receipt result is not an object")
    _require(_canonical(candidate) == _canonical(previous), "permanent receipt result content mismatch")
    # The standalone prior result hash is fixed independently in the manifest;
    # retain this explicit argument so a caller cannot silently substitute one.
    _require(previous_sha == FIXED_SHA256[PREVIOUS_RESULT], "prior result fixed SHA mismatch")


def _extract_output(previous: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(previous.get("status") == "returned", "previous result is not returned")
    _require(previous.get("schema") == "reserved_partial_wholez_v3_result_v1", "previous result schema")
    output = previous.get("output")
    _require(isinstance(output, Mapping), "previous result output missing")
    for key in ("response", "cells", "mu", "joint", "pulse"):
        _require(key in output, "previous result output." + key + " missing")
    fixed_parameters = output.get("fixed_parameters")
    selected_partial_row = output.get("selected_partial_row")
    _require(isinstance(fixed_parameters, Mapping), "previous fixed_parameters missing")
    _require(isinstance(selected_partial_row, Mapping), "previous selected_partial_row missing")
    _require(fixed_parameters.get("P1_consumed") is False, "P1 was consumed")
    _require(selected_partial_row.get("x") == "1.27", "unexpected response point")
    return output


def _sensitivity_mapping(value: Any, label: str) -> dict[str, Any]:
    """Convert the leaf's ordered ``((name, coefficient), ...)`` result."""

    _require(isinstance(value, (list, tuple)), label + " must be an ordered pair tuple")
    pairs: dict[str, Any] = {}
    for index, pair in enumerate(value):
        _require(isinstance(pair, (list, tuple)) and len(pair) == 2, f"{label}[{index}] must be a pair")
        name, coefficient = pair
        _require(isinstance(name, str), f"{label}[{index}] name is not text")
        _require(name not in pairs, f"{label} has duplicate name: {name}")
        pairs[name] = coefficient
    _require(set(pairs) == {"r0_Z", "r1_Z", "r2_Z", "mu_Z"}, label + " names are incomplete")
    return pairs


def _consume(raw: Mapping[str, bytes], hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
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
        _require(isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True, "synthetic fixture failed")
    except BaseException as exc:
        progress["fixture_error"] = _error(exc)
        raise
    progress["fixture_completed_before_input_json_decode"] = True

    progress["input_json_decode_started"] = True
    previous = _load_json_bytes(raw[PREVIOUS_RESULT], PREVIOUS_RESULT)
    receipt = _decode_archive(raw[PERMANENT_RECEIPT], PERMANENT_RECEIPT)
    _receipt_source_result(receipt, previous, hashes[PREVIOUS_RESULT])
    output = _extract_output(previous)

    response_raw = output["response"]
    _require(isinstance(response_raw, list) and len(response_raw) == 3, "output.response must have three intervals")
    response = tuple(_make_interval(c, item, f"output.response[{index}]") for index, item in enumerate(response_raw))
    mu = _make_interval(c, output["mu"], "output.mu")
    joint_raw = output["joint"]
    pulse_raw = output["pulse"]
    _require(isinstance(joint_raw, list) and len(joint_raw) == 2, "output.joint must be [J0,J1]")
    _require(isinstance(pulse_raw, list) and len(pulse_raw) == 2, "output.pulse must be [P0,P1]")
    joint = tuple(_make_interval(c, item, f"output.joint[{index}]") for index, item in enumerate(joint_raw))
    pulse = tuple(_make_interval(c, item, f"output.pulse[{index}]") for index, item in enumerate(pulse_raw))

    cells_raw = output["cells"]
    _require(isinstance(cells_raw, list) and len(cells_raw) == 4, "output.cells must contain four v3 cells")
    evaluate = getattr(leaf, "evaluate_partial_reserved_minus_pulse", None)
    _require(callable(evaluate), "new leaf API is missing")
    unknown_names = ("r0_Z", "r1_Z", "r2_Z", "mu_Z")
    cell_results: list[dict[str, Any]] = []
    actual_leaf_calls = 0
    for index, cell_raw in enumerate(cells_raw):
        _require(isinstance(cell_raw, Mapping), f"output.cells[{index}] is not an object")
        cell_id = cell_raw.get("cell_id", cell_raw.get("source_id", cell_raw.get("id")))
        _require(cell_id is not None, f"output.cells[{index}] has no stable id")
        clipped_Z = _make_interval(c, cell_raw.get("clipped_Z"), f"output.cells[{index}].clipped_Z")
        theta_raw = cell_raw.get("normalized_theta_inputs")
        _require(isinstance(theta_raw, list) and len(theta_raw) >= 2, f"output.cells[{index}] theta orders")
        _require(all(isinstance(order, list) and len(order) == 2 for order in theta_raw[:2]), f"output.cells[{index}] theta shape")
        theta_i = tuple(_make_interval(c, theta_raw[0][component], f"cell[{index}].theta_i[{component}]") for component in range(2))
        theta_i_Z = tuple(_make_interval(c, theta_raw[1][component], f"cell[{index}].theta_i_Z[{component}]") for component in range(2))
        actual_leaf_calls += 1
        progress["actual_leaf_calls"] = actual_leaf_calls
        built = evaluate(
            c,
            response,
            None,
            mu,
            None,
            theta_i,
            theta_i_Z,
            joint[0],
            joint[1],
            pulse[0],
            pulse[1],
            clipped_Z,
        )
        built_fields = _as_fields(built, f"cell[{index}] leaf result")
        for required in ("C0", "Z1_baseline", "Z1", "sensitivity_coefficients", "C0_terms"):
            _require(required in built_fields, f"cell[{index}] leaf result missing {required}")
        _require(built_fields["Z1"] is None, f"cell[{index}] Z1 must be None with unknown derivatives")
        coefficients = _sensitivity_mapping(
            built_fields["sensitivity_coefficients"],
            f"cell[{index}] sensitivity_coefficients",
        )
        coefficient_packets = {
            name: _interval_packet(c, coefficients[name], f"cell[{index}] sensitivity {name}")
            for name in unknown_names
        }
        cell_results.append(
            {
                "source_cell_index": index,
                "cell_id": cell_id,
                "clipped_Z": _interval_packet(c, clipped_Z, f"cell[{index}] clipped_Z"),
                "theta_i": [_interval_packet(c, value, f"cell[{index}] theta_i[{j}]") for j, value in enumerate(theta_i)],
                "theta_i_Z": [_interval_packet(c, value, f"cell[{index}] theta_i_Z[{j}]") for j, value in enumerate(theta_i_Z)],
                "unknown_derivative_names": list(unknown_names),
                "derivative_inputs_omitted_for_baseline": True,
                "derivatives_are_unknown_not_assumed_zero": True,
                "Z1_is_none_with_unknown_derivatives": True,
                "leaf_result": _json_safe(built, context=c, origin=f"cell[{index}] leaf result"),
                "sensitivity_coefficients": coefficient_packets,
            }
        )
    _require(actual_leaf_calls == len(cells_raw), "actual leaf call count mismatch")
    progress["actual_leaf_calls"] = actual_leaf_calls
    return {
        "schema": NAME + "_output_v1",
        "precision_dps": INTERVAL_DPS,
        "previous_result_sha256": hashes[PREVIOUS_RESULT],
        "permanent_receipt_raw_sha256": hashes[PERMANENT_RECEIPT],
        "permanent_receipt_preserves_source_result": True,
        "response": [_interval_packet(c, value, "previous output.response") for value in response],
        "cells": cell_results,
        "mu": _interval_packet(c, mu, "previous output.mu"),
        "joint": [_interval_packet(c, value, "previous output.joint") for value in joint],
        "pulse": [_interval_packet(c, value, "previous output.pulse") for value in pulse],
        "unknown_derivative_names": list(unknown_names),
        "ordinary_derivative_data": {
            "theta_i_Z": True,
            "J_Z": True,
            "P_Z": True,
            "P1_consumed": False,
            "unknown_derivatives_omitted_from_actual_baseline": True,
        },
        "physical_Rp_used": False,
        "function_points_materialized": False,
        "live_owner_claimed": False,
        "actual_leaf_calls": actual_leaf_calls,
    }


def _error(exc: BaseException) -> dict[str, str]:
    return {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}


def main() -> int:
    out: dict[str, Any] = {
        "schema": NAME + "_result_v1",
        "component": NAME,
        "phase": PHASE,
        "status": "failed",
        "actual_leaf_calls": 0,
        "fixture_before_input_json_decode": True,
        "function_points_materialized": False,
        "physical_Rp_used": False,
    }
    result_sha: str | None = None
    try:
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
                "all_raw_inputs_authenticated_before_import_or_decode": True,
                "fixture_before_input_json_decode": True,
                "actual_leaf_calls": 0,
                "function_points_materialized": False,
                "physical_Rp_used": False,
            },
        )
        out.update(
            {
                "source_sha256": source_sha,
                "brief_sha256": brief_sha,
                "entry_sha256": entry_sha,
                "input_hashes": hashes,
                "source_manifest": source,
                "all_raw_inputs_authenticated_before_import_or_decode": True,
            }
        )
        raw = _all_authenticated_bytes(hashes)
        out["all_raw_inputs_authenticated_after_entry"] = True
        output = _consume(raw, hashes, out)
    except BaseException as exc:
        out["error"] = _error(exc)
    else:
        out["status"] = "returned"
        out["output"] = output
    try:
        result_sha = _write_exclusive(RESULT, out)
    except BaseException as exc:
        # Keep the original traceback if result creation itself is unavailable.
        out.setdefault("result_write_error", _error(exc))
    print(
        json.dumps(
            {
                "component": NAME,
                "phase": PHASE,
                "status": out["status"],
                "result_sha256": result_sha,
                "error": out.get("error", {}).get("message") if isinstance(out.get("error"), Mapping) else None,
            },
            sort_keys=True,
        )
    )
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
