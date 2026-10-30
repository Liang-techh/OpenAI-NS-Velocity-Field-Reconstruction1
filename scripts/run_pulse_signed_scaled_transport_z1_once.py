"""One-entry synthetic runner for the caller-owned transport leaf.

The pending brief and source manifest are pre-created inputs.  The source
manifest must authenticate exactly the new leaf and this runner by raw
SHA-256.  At runtime this runner first performs an ``ast.parse`` of the new
leaf, then writes ``entered.json`` with exclusive-create semantics, and only
after that imports and calls the new leaf.  It never imports or calls an old
consumer.  All fixture values are synthetic exact ``Fraction`` values and are
independent hand constants; no source values or physical conclusion are
claimed.
"""

from __future__ import annotations

import ast
import dataclasses
import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "pulse_signed_scaled_transport_z1_v1"
PHASE = COMPONENT + "_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/pulse_signed_scaled_transport_z1.py"
RUNNER_REL = "scripts/run_pulse_signed_scaled_transport_z1_once.py"
LEAF = ROOT / LEAF_REL
RUNNER = ROOT / RUNNER_REL
REQUIRED_INPUTS = frozenset({LEAF_REL, RUNNER_REL})
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


BASE_INPUTS: dict[str, Fraction] = {
    "omega": Fraction(2),
    "omega_Z": Fraction(3),
    "ell_Z": Fraction(1, 3),
    "kappa_Z": Fraction(2, 5),
    "t": Fraction(1, 4),
    "SF": Fraction(3),
    "SF_Z": Fraction(5),
    "dS": Fraction(7),
    "dS_Z": Fraction(11),
    "QyF": Fraction(13),
    "QyF_Z": Fraction(17),
    "dQy": Fraction(19),
    "dQy_Z": Fraction(23),
}

MOVING_INPUTS: dict[str, Fraction] = {
    "ell_xi": Fraction(2),
    "kappa_xi": Fraction(1, 2),
    # The new API requires these source-difference xi jets as part of the
    # all-or-nothing moving block, so keep both explicitly nonzero.
    "SF_xi": Fraction(29),
    "dS_xi": Fraction(31),
    "QyF_xi": Fraction(3),
    "dQy_xi": Fraction(4),
    "xi_endpoint_Z": Fraction(1, 2),
}

HAND_EXPECTED_FIXED = {
    "t2": Fraction(9, 16),
    "B": Fraction(101, 8),
    "B_Z": Fraction(65, 2),
    "H": Fraction(27),
    "H_Z": Fraction(49),
    "density_mantissa": Fraction(101, 8),
    "density_Z_mantissa": Fraction(2873, 48),
    "endpoint_mantissa": Fraction(27),
    "endpoint_Z_mantissa": Fraction(58),
}

HAND_EXPECTED_MOVING = {
    "B_xi": Fraction(643, 8),
    "H_xi": Fraction(103, 4),
    "density_xi_mantissa": Fraction(1047, 8),
    "density_Z_total_mantissa": Fraction(3007, 24),
    "endpoint_xi_mantissa": Fraction(319, 4),
    "endpoint_Z_total_mantissa": Fraction(783, 8),
}

FALSE_FLAGS = {
    "source_ownership": False,
    "physical_scale": False,
    "certified": False,
    "source_ownership_proved": False,
    "physical_scale_proved": False,
    "certification": False,
    "source_identity_proved": False,
    "physical_source_bound": False,
}


class RunnerInputError(RuntimeError):
    """Raised when the pending brief or raw source manifest is not frozen."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any) -> Any:
    """Preserve exact synthetic values in JSON without float coercion."""

    if isinstance(value, Fraction):
        return str(value)
    if dataclasses.is_dataclass(value):
        return {
            field.name: _json_safe(getattr(value, field.name))
            for field in dataclasses.fields(value)
        }
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return {"type": type(value).__name__, "repr": repr(value)}


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _write_once(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_bytes())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _repo_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    candidate = (ROOT / relative).resolve()
    _require(
        candidate.is_relative_to(ROOT.resolve()) and candidate.is_file(),
        label + " is missing or escapes repository: " + normalized,
    )
    return key, candidate


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate exactly the new leaf and this runner by raw bytes."""

    files = source.get("files")
    _require(
        isinstance(files, list) and len(files) == 2,
        "source manifest must contain exactly leaf and runner entries",
    )
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_path(item.get("path"), f"source.files[{index}]")
        _require(key in REQUIRED_INPUTS, "manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "duplicate source manifest path: " + key)
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        hashes[key] = actual
    _require(
        frozenset(hashes) == REQUIRED_INPUTS,
        "manifest must contain exactly leaf and runner raw SHA-256 entries",
    )
    return hashes


def _preflight() -> tuple[dict[str, Any], dict[str, Any], dict[str, str], str, str]:
    _require(DIRECTORY.is_dir(), "component directory is missing")
    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = _load_object(BRIEF, "pending brief")
    source = _load_object(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal guards changed",
    )
    checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    checkpoint_path = checkpoint.get("path")
    _require(
        isinstance(checkpoint_path, str)
        and Path(checkpoint_path).is_absolute(),
        "pending brief source checkpoint path must be absolute",
    )
    _require(
        Path(checkpoint_path).resolve() == SOURCE.resolve(),
        "pending brief source checkpoint path mismatch",
    )
    checkpoint_sha = checkpoint.get("sha256")
    source_sha = _sha256_file(SOURCE)
    _require(
        isinstance(checkpoint_sha, str)
        and SHA256_RE.fullmatch(checkpoint_sha) is not None
        and checkpoint_sha.lower() == source_sha.lower(),
        "pending brief source checkpoint SHA mismatch",
    )
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(
        source.get("schema") == COMPONENT + "_checkpoint_v1",
        "source manifest schema mismatch",
    )
    for key, expected in (
        ("one_entry_only", True),
        ("synthetic_only", True),
        ("no_source_graph_traversal", True),
        ("no_old_function_calls", True),
    ):
        if key in source:
            _require(source[key] is expected, "source manifest guard changed: " + key)
    for key in (
        "source_identity_proved",
        "source_values_computed",
        "physical_scale_proved",
        "certified",
        "historical_error_eliminated",
    ):
        if key in source:
            _require(source[key] is False, "source manifest made forbidden claim: " + key)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, hashlib.sha256(brief_raw).hexdigest(), checkpoint_sha


def _parse_new_leaf() -> None:
    _require(LEAF.is_file(), "new transport leaf is missing")
    try:
        ast.parse(LEAF.read_text(encoding="utf-8"), filename=str(LEAF), mode="exec")
    except (OSError, UnicodeDecodeError, SyntaxError) as exc:
        raise RunnerInputError("new transport leaf failed static ast.parse") from exc


def _import_new_leaf() -> Any:
    spec = importlib.util.spec_from_file_location(
        "pulse_signed_scaled_transport_z1_v1_new_leaf", LEAF
    )
    _require(spec is not None and spec.loader is not None, "new transport leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "leaf result is missing field: " + name)
        return value[name]
    _require(hasattr(value, name), "leaf result is missing field: " + name)
    return getattr(value, name)


def _check_equal(actual: Any, expected: Any, label: str) -> None:
    if isinstance(expected, Fraction):
        _require(isinstance(actual, Fraction), label + " did not retain Fraction type")
    _require(actual == expected, f"{label}: expected {expected}, got {actual}")


def _check_false_flags(value: Any) -> None:
    for name, expected in FALSE_FLAGS.items():
        _check_equal(_field(value, name), expected, "flag " + name)
    _check_equal(_field(value, "status"), "caller_formula_only", "formula status")


def _invoke(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    kwargs: Mapping[str, Any],
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "callable": getattr(function, "__name__", type(function).__name__),
        "invoked": True,
        "input_names": sorted(kwargs),
    }
    calls.append(record)
    try:
        value = function(**dict(kwargs))
    except BaseException as exc:
        record["returned"] = False
        record["passed"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        raise
    record["returned"] = True
    record["passed"] = True
    record["output"] = _json_safe(value)
    return value


def _expect_failure(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    expected_type: type[BaseException],
    kwargs: Mapping[str, Any],
) -> None:
    record: dict[str, Any] = {
        "name": name,
        "callable": getattr(function, "__name__", type(function).__name__),
        "invoked": True,
        "input_names": sorted(kwargs),
        "expected_exception": expected_type.__name__,
    }
    calls.append(record)
    try:
        function(**dict(kwargs))
    except expected_type as exc:
        record["returned"] = False
        record["passed"] = True
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        return
    except BaseException as exc:
        record["returned"] = False
        record["passed"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        raise RunnerInputError(name + " raised the wrong exception") from exc
    record["returned"] = True
    record["passed"] = False
    raise RunnerInputError(name + " unexpectedly returned")


def _result_payload(
    *,
    status: str,
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    entry_sha: str,
    calls: list[dict[str, Any]],
    leaf_imported: bool,
    fixture: Any = None,
    error: BaseException | None = None,
) -> dict[str, Any]:
    flags = {
        "source_identity_proved": False,
        "source_values_computed": False,
        "physical_scale_proved": False,
        "certified": False,
        "old_math_called": False,
        "source_graph_traversed": False,
        "historical_error_eliminated": False,
        "full_goal_complete": False,
    }
    return {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": status,
        "source_checkpoint_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "manifest_paths_exact": sorted(hashes),
        "raw_sha256_authenticated_before_entry": True,
        "leaf_imported": leaf_imported,
        "leaf_invocations": len(calls),
        "new_leaf_only": True,
        "runtime_model_metadata": None,
        "source_manifest": _json_safe(source),
        "fixture": _json_safe(fixture),
        "calls": _json_safe(calls),
        "flags": flags,
        **flags,
        "premises": {
            "synthetic_only": True,
            "fixed_chart_Z_derivatives": True,
            "moving_endpoint_coordinate_supplied": True,
            "source_owned_values": False,
            "physical_scale_claim": False,
            "uncertified": True,
        },
        "error": _json_safe(error),
        "traceback": traceback.format_exc() if error is not None else None,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha, _checkpoint_sha = _preflight()
    _parse_new_leaf()

    # This is intentionally the first state-changing artifact.  No leaf import
    # or call is allowed before this exclusive entry receipt succeeds.
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "manifest_paths_exact": sorted(hashes),
            "leaf_ast_parsed": True,
            "new_leaf_imported": False,
            "leaf_invocations": 0,
            "old_math_called": False,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "certified": False,
            "runtime_model_metadata": None,
        },
    )

    calls: list[dict[str, Any]] = []
    fixture: dict[str, Any] | None = None
    leaf_imported = False
    try:
        module = _import_new_leaf()
        leaf_imported = True
        function = getattr(module, "evaluate_signed_scaled_transport_z1", None)
        _require(callable(function), "new transport leaf callable is missing")
        expected_error = getattr(module, "ScaledTransportInputError", ValueError)
        _require(
            isinstance(expected_error, type) and issubclass(expected_error, BaseException),
            "new transport leaf error type is missing",
        )

        fixed = _invoke(calls, "fixed_hand_constants", function, BASE_INPUTS)
        for field, expected in HAND_EXPECTED_FIXED.items():
            _check_equal(_field(fixed, field), expected, "fixed " + field)
        _check_false_flags(fixed)

        moving_kwargs = dict(BASE_INPUTS)
        moving_kwargs.update(MOVING_INPUTS)
        moving = _invoke(calls, "moving_endpoint_hand_constants", function, moving_kwargs)
        for field, expected in HAND_EXPECTED_FIXED.items():
            _check_equal(_field(moving, field), expected, "moving fixed " + field)
        for field, expected in HAND_EXPECTED_MOVING.items():
            _check_equal(_field(moving, field), expected, "moving " + field)
        _check_false_flags(moving)

        missing_moving = dict(BASE_INPUTS)
        missing_moving["ell_xi"] = MOVING_INPUTS["ell_xi"]
        _expect_failure(
            calls,
            "reject_missing_moving_input",
            function,
            expected_error,
            missing_moving,
        )

        reject_float = dict(BASE_INPUTS)
        reject_float["omega"] = 2.0  # type: ignore[assignment]
        _expect_failure(
            calls,
            "reject_plain_float",
            function,
            expected_error,
            reject_float,
        )

        fixture = {
            "base_inputs": BASE_INPUTS,
            "moving_inputs": MOVING_INPUTS,
            "hand_expected_fixed": HAND_EXPECTED_FIXED,
            "hand_expected_moving": HAND_EXPECTED_MOVING,
            "checks": {
                "fixed_constants": True,
                "moving_endpoint_total": True,
                "missing_moving_rejected": True,
                "plain_float_rejected": True,
                "all_scale_source_certification_flags_false": True,
            },
        }
        result = _result_payload(
            status="returned",
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            brief_sha=brief_sha,
            entry_sha=entry_sha,
            calls=calls,
            leaf_imported=leaf_imported,
            fixture=fixture,
        )
    except BaseException as exc:
        result = _result_payload(
            status="failed",
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            brief_sha=brief_sha,
            entry_sha=entry_sha,
            calls=calls,
            leaf_imported=leaf_imported,
            fixture=fixture,
            error=exc,
        )

    result_sha = _write_once(RESULT, result)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": result["status"],
                "result_sha256": result_sha,
                "leaf_ast_parsed": True,
                "leaf_imported": leaf_imported,
                "leaf_invocations": len(calls),
                "new_leaf_only": True,
                "source_values_computed": False,
                "physical_scale_proved": False,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if result["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
