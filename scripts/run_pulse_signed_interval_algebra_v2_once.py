"""One-entry runner for the exact pulse signed-interval algebra fixture.

The pending source checkpoint contains exactly three raw input hashes: this
runner, the new leaf, and the pinned v1 failure receipt as raw historical
evidence.  The runner authenticates those bytes before loading the leaf, writes
``entered.json`` before the one runner-owned fixture, and writes a complete
result receipt before printing.  It never parses or executes the v1 failure,
imports an older module, or opens a source graph.  The fixture makes several
separately logged synthetic evaluator calls.
"""

from __future__ import annotations

import dataclasses
from fractions import Fraction
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
COMPONENT = "pulse_signed_interval_algebra_v2"
PHASE = COMPONENT + "_pending"
SOURCE_SCHEMA = COMPONENT + "_checkpoint_v1"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "pulse_signed_interval_algebra_v2.py"
RUNNER = ROOT / "scripts" / "run_pulse_signed_interval_algebra_v2_once.py"
LEAF_REL = "src/openai_ns_reconstruction/pulse_signed_interval_algebra_v2.py"
RUNNER_REL = "scripts/run_pulse_signed_interval_algebra_v2_once.py"
PRIOR_FAILURE_REL = "tmp/pulse_signed_interval_algebra_v1/result.json"
PRIOR_FAILURE_SHA256 = "fe27b228b8bea63b4bc19c2c948621bd58825aa9db03473fb65cc6d54dc4c1c1"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

EXPECTED_VALUES = {
    "S_F": Fraction(5, 2),
    "S_F_Z": Fraction(-7, 4),
    "dS": Fraction(3, 1),
    "dS_Z": Fraction(-177, 16),
    "L": Fraction(1, 1),
}
INDEPENDENT_P_VALUES = {
    "S_P": Fraction(11, 2),
    "S_P_Z": Fraction(-205, 16),
}


@dataclasses.dataclass(frozen=True)
class RationalInterval:
    """Small exact rational-endpoint interval used only by the fixture."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        for name in ("lower", "upper"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
                raise TypeError(name + " must be an int or Fraction")
            object.__setattr__(self, name, Fraction(value))
        if self.lower > self.upper:
            raise ValueError("interval lower exceeds upper")

    @classmethod
    def point(cls, value: int | Fraction) -> "RationalInterval":
        value = Fraction(value)
        return cls(value, value)

    @classmethod
    def around(cls, value: int | Fraction, radius: int | Fraction) -> "RationalInterval":
        center = Fraction(value)
        radius = Fraction(radius)
        if radius < 0:
            raise ValueError("interval radius must be nonnegative")
        return cls(center - radius, center + radius)

    @staticmethod
    def _coerce(value: Any) -> "RationalInterval":
        if isinstance(value, RationalInterval):
            return value
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise TypeError("interval arithmetic accepts only exact rational endpoints")
        return RationalInterval.point(value)

    def __add__(self, other: Any) -> "RationalInterval":
        other = self._coerce(other)
        return type(self)(self.lower + other.lower, self.upper + other.upper)

    def __radd__(self, other: Any) -> "RationalInterval":
        return self + other

    def __sub__(self, other: Any) -> "RationalInterval":
        other = self._coerce(other)
        return type(self)(self.lower - other.upper, self.upper - other.lower)

    def __rsub__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) - self

    def __neg__(self) -> "RationalInterval":
        return type(self)(-self.upper, -self.lower)

    def __mul__(self, other: Any) -> "RationalInterval":
        other = self._coerce(other)
        values = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return type(self)(min(values), max(values))

    def __rmul__(self, other: Any) -> "RationalInterval":
        return self * other

    def __truediv__(self, other: Any) -> "RationalInterval":
        other = self._coerce(other)
        if other.lower <= 0 <= other.upper:
            raise ZeroDivisionError("interval denominator contains zero")
        reciprocal_values = (
            Fraction(1, 1) / other.lower,
            Fraction(1, 1) / other.upper,
        )
        reciprocal = type(self)(min(reciprocal_values), max(reciprocal_values))
        return self * reciprocal

    def __rtruediv__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) / self

    def contains(self, value: int | Fraction) -> bool:
        value = Fraction(value)
        return self.lower <= value <= self.upper

    def __str__(self) -> str:
        return f"[{self.lower},{self.upper}]"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _json_safe(value: Any) -> Any:
    """Keep exact fixture rationals and interval endpoints JSON-visible."""

    if isinstance(value, Fraction):
        return {
            "type": "Fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
        }
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _json_safe(dataclasses.asdict(value))
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return {"type": type(value).__name__, "fields": _json_safe(fields)}
    return str(value)


def _write_exclusive(path: Path, value: Any) -> str:
    raw = _json_bytes(_json_safe(value))
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _relative_path(raw_path: Any, label: str) -> str:
    _require(isinstance(raw_path, str) and raw_path, label + " path missing")
    candidate = Path(raw_path)
    if candidate.is_absolute():
        try:
            candidate = candidate.resolve().relative_to(ROOT.resolve())
        except ValueError as exc:
            raise RuntimeError(label + " path is outside repository") from exc
    _require(".." not in candidate.parts, label + " path escapes repository")
    return candidate.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Read direct path/hash records without traversing any source artifact."""

    records: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    _require(isinstance(files, list) and files, "source files manifest missing")
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        records.append((item.get("path"), item.get("sha256"), f"source.files[{index}]"))
    for key in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256", "code_sha256"):
        table = source.get(key)
        if table is None:
            continue
        _require(isinstance(table, Mapping), "source." + key + " is not a mapping")
        records.extend((name, digest, f"source.{key}[{name!r}]") for name, digest in table.items())

    result: dict[str, str] = {}
    for raw_path, raw_digest, label in records:
        relative = _relative_path(raw_path, label)
        _require(isinstance(raw_digest, str) and SHA256_RE.fullmatch(raw_digest) is not None,
                 label + " has no valid SHA-256")
        digest = raw_digest.lower()
        _require(relative not in result or result[relative] == digest,
                 "conflicting source SHA-256: " + relative)
        result[relative] = digest
    _require(set(result) == {LEAF_REL, RUNNER_REL, PRIOR_FAILURE_REL},
             "source must hash exactly v2 leaf, v2 runner, and v1 failure")
    _require(result[PRIOR_FAILURE_REL] == PRIOR_FAILURE_SHA256,
             "v1 failure SHA pin changed")
    return result


def _brief_checkpoint(brief: Mapping[str, Any]) -> tuple[str, str]:
    for key in (f"{COMPONENT}_source_checkpoint", "source_checkpoint"):
        value = brief.get(key)
        if isinstance(value, Mapping):
            path, digest = value.get("path"), value.get("sha256")
            if isinstance(path, str) and isinstance(digest, str):
                return path, digest
    path = brief.get("source_checkpoint_path")
    digest = brief.get("source_checkpoint_sha256")
    _require(isinstance(path, str) and isinstance(digest, str),
             "pending source checkpoint pointer missing")
    return path, digest


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTERED.exists() and not RESULT.exists(), "already entered; never re-enter")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source checkpoint is missing")

    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    try:
        brief = json.loads(brief_raw.decode("utf-8"))
        source = json.loads(source_raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError("brief/source checkpoint is not valid JSON") from exc
    _require(isinstance(brief, Mapping), "pending brief is not an object")
    _require(isinstance(source, Mapping), "source checkpoint is not an object")

    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    for key, expected in (
        ("goal_active", True),
        ("full_goal_complete", False),
        ("do_not_reenter", True),
    ):
        if key in brief:
            _require(brief.get(key) is expected, "pending brief guard failed: " + key)

    _require(source.get("schema") == SOURCE_SCHEMA, "source checkpoint schema mismatch")
    _require(source.get("component") == COMPONENT, "source checkpoint component mismatch")
    _require(source.get("one_entry_only") is True, "source checkpoint is not one-entry-only")
    _require(source.get("no_source_graph_traversal") is True,
             "source graph traversal guard missing")

    checkpoint_path, checkpoint_expected = _brief_checkpoint(brief)
    checkpoint_relative = _relative_path(checkpoint_path, "brief source checkpoint")
    _require((ROOT / checkpoint_relative).resolve() == SOURCE.resolve(),
             "brief points to another source checkpoint")
    _require(SHA256_RE.fullmatch(checkpoint_expected) is not None,
             "brief checkpoint SHA malformed")
    source_sha = hashlib.sha256(source_raw).hexdigest()
    _require(source_sha == checkpoint_expected.lower(), "source checkpoint SHA mismatch")

    input_hashes = _manifest_hashes(source)
    # Every listed byte is authenticated before the leaf can be loaded.
    for relative, expected in input_hashes.items():
        path = ROOT / relative
        _require(path.is_file(), "frozen input is missing: " + relative)
        _require(sha256(path) == expected, "frozen raw input changed: " + relative)
    return dict(source), input_hashes, source_sha, hashlib.sha256(brief_raw).hexdigest()


def _fraction(value: Any) -> Fraction | None:
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value, 1)
    if isinstance(value, str):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError):
            return None
    if isinstance(value, Mapping):
        numerator = value.get("numerator")
        denominator = value.get("denominator")
        if type(numerator) is int and type(denominator) is int and denominator != 0:
            try:
                return Fraction(numerator, denominator)
            except ZeroDivisionError:
                return None
    return None


def _contains_builtin_float(value: Any) -> bool:
    if type(value) is float:
        return True
    if isinstance(value, Mapping):
        return any(_contains_builtin_float(item) for item in value.values())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_builtin_float(item) for item in value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _contains_builtin_float(dataclasses.asdict(value))
    fields = getattr(value, "__dict__", None)
    return isinstance(fields, dict) and _contains_builtin_float(fields)


def _result_fields(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: getattr(value, field.name)
            for field in dataclasses.fields(value)
        }
    fields = getattr(value, "__dict__", None)
    _require(isinstance(fields, dict), "leaf result has no named fields")
    return fields


def _check_exact_result(value: Any, *, require_fraction: bool) -> dict[str, Any]:
    fields = _result_fields(value)
    result: dict[str, Any] = {}
    for name, expected in EXPECTED_VALUES.items():
        observed_value = fields.get(name)
        observed = _fraction(observed_value)
        _require(observed == expected, "exact value mismatch: " + name)
        if require_fraction:
            _require(type(observed_value) is Fraction, "result is not a Fraction: " + name)
        result[name] = observed_value
    return result


def _invoke_leaf(
    evaluate: Any,
    calls: list[dict[str, Any]],
    *,
    name: str,
    frame: Mapping[str, Any],
    signed: Mapping[str, Any],
    Z: Any,
    delta: Any,
    delta_Z: Any,
    gate: Any,
    expected_exception: type[BaseException] | None = None,
    require_fraction: bool = False,
    check_exact: bool = True,
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "leaf_function": "evaluate_signed_interval_algebra",
        "invoked": False,
        "expected_exception": expected_exception.__name__ if expected_exception else None,
    }
    calls.append(record)
    try:
        record["invoked"] = True
        result = evaluate(
            frame,
            signed,
            Z=Z,
            delta=delta,
            delta_Z=delta_Z,
            denominator_nonzero=gate,
        )
        record["returned"] = True
        record["output"] = _json_safe(result)
        if expected_exception is not None:
            record["passed"] = False
            raise RuntimeError(name + " unexpectedly returned")
        if check_exact:
            fields = _check_exact_result(result, require_fraction=require_fraction)
            record["exact_values"] = _json_safe(fields)
        record["passed"] = True
        return result
    except BaseException as exc:
        record["returned"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        if expected_exception is None:
            record["passed"] = False
            raise
        record["passed"] = isinstance(exc, expected_exception)
        if not record["passed"]:
            raise
        return None


def _run_synthetic_fixture(leaf: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    """Call the new evaluator only through this one runner-owned fixture."""

    evaluate = getattr(leaf, "evaluate_signed_interval_algebra", None)
    _require(callable(evaluate), "new leaf evaluator is missing")
    frame_names = ("U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ")
    signed_names = ("dU", "dU_Z", "dQ", "dQ_Z", "dQ_ZZ", "dQ_y", "dQ_yZ")
    fraction_frame = {name: Fraction(1, 1) for name in frame_names}
    fraction_signed = {name: Fraction(1, 1) for name in signed_names}
    integer_frame = {name: 1 for name in frame_names}
    integer_signed = {name: 1 for name in signed_names}
    parameters = {
        "Z": Fraction(1, 2),
        "delta": Fraction(0, 1),
        "delta_Z": Fraction(1, 2),
    }

    fraction_result = _invoke_leaf(
        evaluate,
        calls,
        name="fraction_exact",
        frame=fraction_frame,
        signed=fraction_signed,
        **parameters,
        gate=lambda _L: True,
        require_fraction=True,
    )
    integer_result = _invoke_leaf(
        evaluate,
        calls,
        name="builtin_int_exact",
        frame=integer_frame,
        signed=integer_signed,
        **parameters,
        gate=lambda _L: True,
        require_fraction=True,
    )
    zero_signed = {name: Fraction(0, 1) for name in signed_names}
    zero_delta_result = _invoke_leaf(
        evaluate,
        calls,
        name="direct_signed_delta_zero",
        frame=fraction_frame,
        signed=zero_signed,
        **parameters,
        gate=lambda _L: True,
        check_exact=False,
    )
    zero_delta_fields = _result_fields(zero_delta_result)
    zero_delta_values = {
        "S_F": EXPECTED_VALUES["S_F"],
        "S_F_Z": EXPECTED_VALUES["S_F_Z"],
        "dS": Fraction(0, 1),
        "dS_Z": Fraction(0, 1),
        "L": EXPECTED_VALUES["L"],
    }
    zero_delta_checks = {
        name: _fraction(zero_delta_fields.get(name)) == expected
        and type(zero_delta_fields.get(name)) is Fraction
        for name, expected in zero_delta_values.items()
    }
    _require(all(zero_delta_checks.values()), "zero signed delta result mismatch")
    _invoke_leaf(
        evaluate,
        calls,
        name="denominator_gate_false",
        frame=fraction_frame,
        signed=fraction_signed,
        **parameters,
        gate=lambda _L: False,
        expected_exception=ValueError,
    )
    _invoke_leaf(
        evaluate,
        calls,
        name="denominator_gate_nonbool",
        frame=fraction_frame,
        signed=fraction_signed,
        **parameters,
        gate=lambda _L: 1,
        expected_exception=ValueError,
    )
    float_frame = dict(fraction_frame)
    float_frame["U"] = 1.0
    _invoke_leaf(
        evaluate,
        calls,
        name="builtin_float_rejected",
        frame=float_frame,
        signed=fraction_signed,
        **parameters,
        gate=lambda _L: True,
        expected_exception=TypeError,
    )
    _invoke_leaf(
        evaluate,
        calls,
        name="zero_L_gate_rejected",
        frame=fraction_frame,
        signed=fraction_signed,
        Z=Fraction(1, 1),
        delta=Fraction(1, 1),
        delta_Z=Fraction(0, 1),
        gate=lambda L: L != 0,
        expected_exception=ValueError,
    )

    interval_frame = {
        name: RationalInterval.around(1, Fraction(1, 100)) for name in frame_names
    }
    interval_signed = {
        name: RationalInterval.around(1, Fraction(1, 100)) for name in signed_names
    }
    interval_result = _invoke_leaf(
        evaluate,
        calls,
        name="rational_interval_backend",
        frame=interval_frame,
        signed=interval_signed,
        Z=RationalInterval.around(Fraction(1, 2), Fraction(1, 100)),
        delta=RationalInterval.around(0, Fraction(1, 100)),
        delta_Z=RationalInterval.around(Fraction(1, 2), Fraction(1, 100)),
        gate=lambda _L: True,
        check_exact=False,
    )
    interval_fields = _result_fields(interval_result)
    interval_contains = {
        name: isinstance(interval_fields.get(name), RationalInterval)
        and interval_fields[name].contains(expected)
        for name, expected in EXPECTED_VALUES.items()
    }
    _require(all(interval_contains.values()), "interval does not contain independent point")
    _invoke_leaf(
        evaluate,
        calls,
        name="interval_L_contains_zero_rejected",
        frame=interval_frame,
        signed=interval_signed,
        Z=RationalInterval.point(Fraction(1, 2)),
        delta=RationalInterval(Fraction(3, 1), Fraction(5, 1)),
        delta_Z=RationalInterval.point(Fraction(1, 2)),
        gate=lambda L: L.lower > 0 or L.upper < 0,
        expected_exception=ValueError,
    )

    fraction_fields = _result_fields(fraction_result)
    independent_differences = {
        "dS": INDEPENDENT_P_VALUES["S_P"] - EXPECTED_VALUES["S_F"],
        "dS_Z": INDEPENDENT_P_VALUES["S_P_Z"] - EXPECTED_VALUES["S_F_Z"],
    }
    hand_calculation_checks = {
        "S_P_matches_independent_hand_value": INDEPENDENT_P_VALUES["S_P"] == Fraction(11, 2),
        "S_P_Z_matches_independent_hand_value": INDEPENDENT_P_VALUES["S_P_Z"] == Fraction(-205, 16),
        "hand_difference_dS_is_three": independent_differences["dS"] == Fraction(3, 1),
        "hand_difference_dS_Z_is_minus_177_over_16": independent_differences["dS_Z"] == Fraction(-177, 16),
        "signed_dS_equals_hand_difference": _fraction(fraction_fields.get("dS")) == independent_differences["dS"],
        "signed_dS_Z_equals_hand_difference": _fraction(fraction_fields.get("dS_Z")) == independent_differences["dS_Z"],
    }
    _require(all(hand_calculation_checks.values()), "independent P hand calculation mismatch")

    checks = {
        "fraction_values_match_hand_calculation": True,
        "fraction_outputs_are_exact": all(
            type(_result_fields(fraction_result).get(name)) is Fraction
            for name in EXPECTED_VALUES
        ),
        "builtin_int_values_match_hand_calculation": True,
        "builtin_int_outputs_are_exact": all(
            type(_result_fields(integer_result).get(name)) is Fraction
            for name in EXPECTED_VALUES
        ),
        "direct_signed_delta_zero": all(zero_delta_checks.values()),
        "denominator_gate_false_rejected": calls[3]["passed"],
        "denominator_gate_nonbool_rejected": calls[4]["passed"],
        "builtin_float_rejected": calls[5]["passed"],
        "zero_L_gate_rejected": calls[6]["passed"],
        "interval_L_contains_zero_rejected": calls[8]["passed"],
        **hand_calculation_checks,
        "interval_has_rational_endpoints": all(
            isinstance(interval_fields[name], RationalInterval)
            and type(interval_fields[name].lower) is Fraction
            and type(interval_fields[name].upper) is Fraction
            for name in EXPECTED_VALUES
        ),
        "interval_contains_independent_point": all(interval_contains.values()),
        "no_builtin_float_exposed": not _contains_builtin_float(
            {"fraction": fraction_result, "integer": integer_result, "interval": interval_result}
        ),
    }
    passed = all(checks.values()) and len(calls) == 9 and all(item["invoked"] for item in calls)
    return {
        "schema": COMPONENT + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_checks_passed" if passed else "synthetic_checks_failed",
        "passed": passed,
        "checks_passed": passed,
        "synthetic_only": True,
        "checks": checks,
        "call_count": len(calls),
        "calls": _json_safe(calls),
        "fraction_result": _json_safe(fraction_result),
        "integer_result": _json_safe(integer_result),
        "zero_delta_result": _json_safe(zero_delta_result),
        "interval_result": _json_safe(interval_result),
        "interval_contains_independent_point": interval_contains,
        "expected_values": {name: str(value) for name, value in EXPECTED_VALUES.items()},
        "independent_P_values": {name: str(value) for name, value in INDEPENDENT_P_VALUES.items()},
        "independent_signed_differences": {
            name: str(value) for name, value in independent_differences.items()
        },
        "zero_delta_checks": zero_delta_checks,
        "hand_calculation_checks": hand_calculation_checks,
        "rounding_certified": False,
        "certified": False,
        "source_identity_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_value_authentication": False,
        "source_values_authenticated": False,
        "source_graph_traversed": False,
        "actual_source_values": False,
    }


def _validate_fixture(value: Any) -> None:
    _require(isinstance(value, Mapping), "runner fixture did not return an object")
    _require(value.get("passed") is True and value.get("checks_passed") is True,
             "runner fixture checks did not pass")
    _require(value.get("synthetic_only") is True, "runner fixture is not synthetic-only")
    _require(value.get("call_count") == 9, "fixture call count changed")
    calls = value.get("calls")
    _require(isinstance(calls, list) and len(calls) == 9, "fixture call log incomplete")
    _require(all(isinstance(item, Mapping) and item.get("invoked") is True for item in calls),
             "fixture call log has an uninvoked call")
    _require(not _contains_builtin_float(value), "fixture exposed a built-in float")
    for key in (
        "certified",
        "rounding_certified",
        "source_identity_proved",
        "source_values_selected",
        "source_values_computed",
        "source_value_authentication",
        "source_values_authenticated",
    ):
        if key in value:
            _require(value[key] is False, "fixture made forbidden claim: " + key)
    checks = value.get("checks")
    _require(isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()),
             "fixture contains a failed check")


def main() -> int:
    _, input_hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTERED,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": dict(input_hashes),
            "fixture_invocations": 0,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
        },
    )

    fixture_calls = 0
    leaf_imported = False
    fixture_audit: list[dict[str, Any]] = []
    try:
        spec = importlib.util.spec_from_file_location(
            "pulse_signed_interval_algebra_v2_new_leaf", LEAF
        )
        _require(spec is not None and spec.loader is not None, "new leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        leaf_imported = True

        fixture_calls = 1
        value = _run_synthetic_fixture(leaf, fixture_audit)
        _validate_fixture(value)
        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": dict(input_hashes),
            "leaf_imported": leaf_imported,
            "fixture_invocations": fixture_calls,
            "consumer_calls": fixture_calls,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "source_value_authentication": False,
            "source_values_authenticated": False,
            "rounding_certified": False,
            "source_identity_proved": False,
            "certified": False,
            "source_receipt": {
                "raw_input_sha_authentication_only": True,
                "source_values_selected": False,
                "source_values_computed": False,
                "source_value_authentication": False,
                "source_identity_proved": False,
            },
            "runtime_model_metadata": None,
            "expected_fixture_contract": {
                "Z": "1/2",
                "delta": "0",
                "delta_Z": "1/2",
                "fourteen_jets": "1",
                "S_F": "5/2",
                "S_F_Z": "-7/4",
                "dS": "3",
                "dS_Z": "-177/16",
                "L": "1",
                "S_P": "11/2",
                "S_P_Z": "-205/16",
                "direct_signed_delta_zero": {"dS": "0", "dS_Z": "0"},
                "required_rejections": [
                    "denominator gate false",
                    "denominator gate non-bool",
                    "built-in float",
                    "zero L",
                    "L interval containing zero",
                ],
                "optional_interval_backend": "exact rational endpoints and independent point containment only",
                "rounding_certification": False,
            },
            "error": None,
            "result": _json_safe(value),
            "fixture_call_log": _json_safe(fixture_audit),
        }
    except BaseException as exc:  # preserve a complete failure artifact
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": dict(input_hashes),
            "leaf_imported": leaf_imported,
            "fixture_invocations": fixture_calls,
            "consumer_calls": fixture_calls,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "source_value_authentication": False,
            "source_values_authenticated": False,
            "rounding_certified": False,
            "source_identity_proved": False,
            "certified": False,
            "source_receipt": {
                "raw_input_sha_authentication_only": True,
                "source_values_selected": False,
                "source_values_computed": False,
                "source_value_authentication": False,
                "source_identity_proved": False,
            },
            "runtime_model_metadata": None,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
            "result": None,
            "fixture_call_log": _json_safe(fixture_audit),
        }

    result_sha = _write_exclusive(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "source_checkpoint_sha256": source_sha,
                "entry_sha256": entry_sha,
                "input_hashes": input_hashes,
                "fixture_invocations": fixture_calls,
                "rounding_certified": False,
                "source_identity_proved": False,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
