"""Enter the synthetic F-gap native-radial-Q jet fixture exactly once.

The coordinated leaf interface is
``reconstruct_f_gap_fixed_xi_q_jets`` with positional ordinary-Z rows in the
following order::

    U0, U1, U2,
    m0, m1, m2, m3,
    Uy0, Uy1,
    my0, my1, my2,
    delta0, delta1, delta2, Z,
    *, denominator_nonzero_gate=Callable[[L], bool], delta_y_zero_premise=True

Here ``m1``, ``m2`` and ``m3`` are the ordinary derivatives
``m_Z``, ``m_ZZ`` and ``m_ZZZ``; they are not Taylor coefficients.  The leaf
returns seven fields in the fixed-``xi`` order
``(U, U_Z, Q, Q_Z, Q_ZZ, Q_y, Q_yZ)``.

This runner owns only an independent exact-rational hand fixture and a small
rational-endpoint interval fixture.  It does not read a source archive, import
an ancestor, authenticate source values, or claim that the fixture values are
real F values.  A pending manifest must contain exactly the new leaf and this
runner, each pinned by its raw SHA-256.  The entry receipt is written before
the leaf is imported; after entry, both successful and failed fixture results
are persisted as a complete ``result.json`` receipt.
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
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_q_jets_from_source_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "f_gap_q_jets_from_source.py"
RUNNER = ROOT / "scripts" / "run_f_gap_q_jets_from_source_once.py"
LEAF_REL = "src/openai_ns_reconstruction/f_gap_q_jets_from_source.py"
RUNNER_REL = "scripts/run_f_gap_q_jets_from_source_once.py"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


class RationalInterval:
    """Closed interval with exact rational endpoints for the fixture only."""

    __slots__ = ("lower", "upper")

    def __init__(self, lower: int | Fraction, upper: int | Fraction) -> None:
        if isinstance(lower, bool) or isinstance(upper, bool):
            raise TypeError("interval endpoints must be exact rationals")
        if not isinstance(lower, (int, Fraction)) or not isinstance(upper, (int, Fraction)):
            raise TypeError("interval endpoints must be exact rationals")
        self.lower = Fraction(lower)
        self.upper = Fraction(upper)
        if self.lower > self.upper:
            raise ValueError("interval lower exceeds upper")

    @classmethod
    def around(cls, value: int | Fraction, radius: int | Fraction) -> "RationalInterval":
        value = Fraction(value)
        radius = Fraction(radius)
        if radius < 0:
            raise ValueError("interval radius must be nonnegative")
        return cls(value - radius, value + radius)

    @classmethod
    def point(cls, value: int | Fraction) -> "RationalInterval":
        value = Fraction(value)
        return cls(value, value)

    @staticmethod
    def _coerce(value: Any) -> "RationalInterval":
        if isinstance(value, RationalInterval):
            return value
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise TypeError("interval arithmetic accepts only exact rationals")
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
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return type(self)(min(products), max(products))

    def __rmul__(self, other: Any) -> "RationalInterval":
        return self * other

    def __truediv__(self, other: Any) -> "RationalInterval":
        other = self._coerce(other)
        if other.lower <= 0 <= other.upper:
            raise ZeroDivisionError("interval denominator contains zero")
        reciprocal = (Fraction(1, 1) / other.lower, Fraction(1, 1) / other.upper)
        return self * type(self)(min(reciprocal), max(reciprocal))

    def __rtruediv__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) / self

    def contains(self, value: int | Fraction) -> bool:
        value = Fraction(value)
        return self.lower <= value <= self.upper


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any) -> Any:
    if isinstance(value, Fraction):
        return {
            "type": "Fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
        }
    if isinstance(value, RationalInterval):
        return {
            "type": "RationalInterval",
            "lower": _json_safe(value.lower),
            "upper": _json_safe(value.upper),
        }
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _json_safe(
            {field.name: getattr(value, field.name) for field in dataclasses.fields(value)}
        )
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return _json_safe(fields)
    if hasattr(value, "lower") and hasattr(value, "upper"):
        return {"lower": _json_safe(value.lower), "upper": _json_safe(value.upper)}
    return repr(value)


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


def _write_exclusive(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
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
            raise RunnerInputError(label + " path is outside repository") from exc
    _require(".." not in candidate.parts, label + " path escapes repository")
    return candidate.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list), "source files manifest missing")
    _require(len(files) == 2, "source manifest must contain exactly leaf plus runner")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        relative = _relative_path(item.get("path"), f"source.files[{index}]")
        _require(
            relative in {LEAF_REL, RUNNER_REL},
            "manifest contains an out-of-scope file: " + relative,
        )
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            f"source.files[{index}] has no valid SHA-256",
        )
        _require(relative not in hashes, "duplicate source manifest path: " + relative)
        path = ROOT / relative
        _require(path.is_file(), "frozen input is missing: " + relative)
        actual = _sha(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + relative)
        hashes[relative] = actual
    _require(set(hashes) == {LEAF_REL, RUNNER_REL}, "manifest must be exactly leaf plus runner")
    return hashes


def _brief_checkpoint(brief: Mapping[str, Any]) -> tuple[str, str]:
    for key in (f"{COMPONENT}_source_checkpoint", "source_checkpoint"):
        value = brief.get(key)
        if isinstance(value, Mapping):
            path, digest = value.get("path"), value.get("sha256")
            if isinstance(path, str) and isinstance(digest, str):
                return path, digest
    path = brief.get("source_checkpoint_path")
    digest = brief.get("source_checkpoint_sha256")
    _require(
        isinstance(path, str) and isinstance(digest, str),
        "pending source checkpoint pointer missing",
    )
    return path, digest


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    """Validate the frozen brief/manifest before any new-leaf import."""

    _require(not ENTERED.exists() and not RESULT.exists(), "already entered; never re-enter")
    _require(DIRECTORY.is_dir(), "frozen component directory is missing")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    try:
        brief = json.loads(brief_raw.decode("utf-8"))
        source = json.loads(source_raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError("brief/source manifest is not valid JSON") from exc
    _require(isinstance(brief, Mapping), "pending brief is not an object")
    _require(isinstance(source, Mapping), "source manifest is not an object")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief guard failed",
    )
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("one_entry_only") is True, "one-entry guard missing")
    _require(source.get("synthetic_only") is True, "synthetic-only guard missing")
    _require(
        source.get("no_source_graph_traversal") is True,
        "no-source-graph-traversal guard missing",
    )
    _require(
        source.get("no_ancestor_replay") is True
        or source.get("no_old_ancestor_replay") is True,
        "no-ancestor-replay guard missing",
    )
    checkpoint_path, checkpoint_expected = _brief_checkpoint(brief)
    checkpoint_relative = _relative_path(checkpoint_path, "brief source checkpoint")
    _require((ROOT / checkpoint_relative).resolve() == SOURCE.resolve(), "brief points elsewhere")
    _require(
        SHA256_RE.fullmatch(checkpoint_expected) is not None,
        "source checkpoint SHA malformed",
    )
    source_sha = hashlib.sha256(source_raw).hexdigest()
    _require(source_sha == checkpoint_expected.lower(), "source checkpoint SHA mismatch")
    hashes = _manifest_hashes(source)
    _require(LEAF.is_file(), "new Q-jet leaf is missing")

    # This component has no source-value input.  If a caller adds any of these
    # semantic fields, a positive claim is forbidden rather than interpreted.
    for key in (
        "source_identity_proved",
        "source_values_selected",
        "source_values_computed",
        "source_value_authentication",
        "authentication_proved",
        "numeric_use_admissible",
        "certified",
    ):
        if key in source:
            _require(source[key] is False, "manifest made forbidden source claim: " + key)
    return dict(source), hashes, source_sha, hashlib.sha256(brief_raw).hexdigest()


def _import_leaf() -> Any:
    spec = importlib.util.spec_from_file_location("f_gap_q_jets_from_source_v1_new_leaf", LEAF)
    _require(spec is not None and spec.loader is not None, "new Q-jet leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _result_fields(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: getattr(value, field.name)
            for field in dataclasses.fields(value)
        }
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return fields
    raise RunnerInputError("leaf result has no named fields")


def _result_tuple(value: Any) -> tuple[Any, ...]:
    as_tuple = getattr(value, "as_tuple", None)
    if callable(as_tuple):
        result = tuple(as_tuple())
    else:
        fields = _result_fields(value)
        result = tuple(fields[name] for name in (
            "U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ"
        ))
    _require(len(result) == 7, "leaf did not return seven Q jets")
    return result


def _as_interval(value: Any) -> RationalInterval | None:
    if isinstance(value, RationalInterval):
        return value
    fields = _result_fields(value) if value is not None else None
    if fields is None:
        return None
    lower = fields.get("lower")
    upper = fields.get("upper")
    if isinstance(lower, (int, Fraction)) and not isinstance(lower, bool) and isinstance(upper, (int, Fraction)) and not isinstance(upper, bool):
        return RationalInterval(lower, upper)
    return None


def _contains_builtin_float(value: Any) -> bool:
    if type(value) is float:
        return True
    if isinstance(value, Mapping):
        return any(_contains_builtin_float(item) for item in value.values())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_builtin_float(item) for item in value)
    if isinstance(value, RationalInterval):
        return _contains_builtin_float(value.lower) or _contains_builtin_float(value.upper)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(
            _contains_builtin_float(getattr(value, field.name))
            for field in dataclasses.fields(value)
        )
    fields = getattr(value, "__dict__", None)
    return isinstance(fields, dict) and _contains_builtin_float(fields)


def _contains_midpoint_field(value: Any) -> bool:
    if isinstance(value, Mapping):
        return any(
            "midpoint" in str(key).lower() or _contains_midpoint_field(item)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_midpoint_field(item) for item in value)
    if isinstance(value, RationalInterval):
        return False
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(
            "midpoint" in field.name.lower()
            or _contains_midpoint_field(getattr(value, field.name))
            for field in dataclasses.fields(value)
        )
    fields = getattr(value, "__dict__", None)
    return isinstance(fields, dict) and _contains_midpoint_field(fields)


def _nonzero_gate(value: Any) -> bool:
    """Return a built-in bool only when the exact/interval denominator excludes zero."""

    if isinstance(value, RationalInterval):
        return bool(value.lower > 0 or value.upper < 0)
    if isinstance(value, Fraction):
        return bool(value != 0)
    return False


def _canonical_inputs() -> tuple[Fraction, ...]:
    # Deliberately positional: m1/m2/m3 must reach the leaf as m_Z/m_ZZ/m_ZZZ.
    return (
        Fraction(1), Fraction(2), Fraction(3),
        Fraction(2), Fraction(1), Fraction(3), Fraction(4),
        Fraction(5), Fraction(6),
        Fraction(7), Fraction(8), Fraction(9),
        Fraction(1, 2), Fraction(1, 3), Fraction(1, 4), Fraction(1, 2),
    )


HAND_EXPECTED = (
    Fraction(1),
    Fraction(2),
    Fraction(-2, 7),
    Fraction(40, 21),
    Fraction(9019, 441),
    Fraction(-22, 7),
    Fraction(266, 21),
)

HAND_CALCULATION = {
    "L": "7/8",
    "L_Z": "-7/12",
    "L_ZZ": "-83/48",
    "N": "-1/4",
    "N_Z": "11/6",
    "N_ZZ": "97/6",
    "N_y": "-11/4",
    "N_yZ": "155/12",
    "returned_order": ["U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ"],
    "jets": ["1", "2", "-2/7", "40/21", "9019/441", "-22/7", "266/21"],
}


def _record_success(
    evaluator: Callable[..., Any],
    calls: list[dict[str, Any]],
    *,
    name: str,
    values: tuple[Any, ...],
    expected: tuple[Fraction, ...] | None = None,
    interval: bool = False,
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "invoked": False,
        "input_scope": "synthetic_exact_rational" if not interval else "synthetic_rational_interval",
        "input_order": [
            "U0", "U1", "U2", "m0", "m1", "m2", "m3",
            "Uy0", "Uy1", "my0", "my1", "my2", "delta0", "delta1", "delta2", "Z",
        ],
    }
    calls.append(record)
    try:
        record["invoked"] = True
        output = evaluator(
            *values,
            denominator_nonzero_gate=_nonzero_gate,
            delta_y_zero_premise=True,
        )
        result = _result_tuple(output)
        record["returned"] = True
        record["output"] = _json_safe(output)
        _require(not _contains_builtin_float(output), name + " exposed a built-in float")
        _require(not _contains_midpoint_field(output), name + " exposed an interval midpoint")
        if interval:
            _require(expected is not None, name + " interval expected values missing")
            for observed, expected_value in zip(result, expected):
                enclosed = _as_interval(observed)
                _require(enclosed is not None, name + " did not return rational intervals")
                _require(enclosed.contains(expected_value), name + " interval misses hand value")
            record["interval_contains_hand_value"] = True
        else:
            _require(expected is not None, name + " exact expected values missing")
            _require(all(type(value) is Fraction for value in result), name + " is not exact Fraction output")
            _require(result == expected, name + " differs from independent hand calculation")
            record["exact_matches_hand_calculation"] = True
        record["passed"] = True
        return output
    except BaseException as exc:
        record["returned"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["passed"] = False
        raise


def _record_expected_failure(
    evaluator: Callable[..., Any],
    calls: list[dict[str, Any]],
    *,
    name: str,
    values: tuple[Any, ...],
    gate: Any,
    delta_y_zero_premise: bool = True,
    expected_exception: tuple[type[BaseException], ...],
) -> None:
    record: dict[str, Any] = {
        "name": name,
        "invoked": False,
        "expected_exception": [item.__name__ for item in expected_exception],
        "input_scope": "negative_synthetic_guard",
    }
    calls.append(record)
    try:
        record["invoked"] = True
        evaluator(
            *values,
            denominator_nonzero_gate=gate,
            delta_y_zero_premise=delta_y_zero_premise,
        )
    except BaseException as exc:
        record["returned"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["passed"] = isinstance(exc, expected_exception)
        if not record["passed"]:
            raise
        return
    raise RunnerInputError(name + " unexpectedly returned")


def _run_fixture(module: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    evaluator = getattr(module, "reconstruct_f_gap_fixed_xi_q_jets", None)
    _require(callable(evaluator), "new leaf callable is missing: reconstruct_f_gap_fixed_xi_q_jets")

    exact_inputs = _canonical_inputs()
    exact = _record_success(
        evaluator,
        calls,
        name="exact_fraction_hand_fixture",
        values=exact_inputs,
        expected=HAND_EXPECTED,
    )

    interval_inputs = tuple(
        RationalInterval.around(value, Fraction(1, 1000)) for value in exact_inputs
    )
    interval = _record_success(
        evaluator,
        calls,
        name="rational_endpoint_interval_fixture",
        values=interval_inputs,
        expected=HAND_EXPECTED,
        interval=True,
    )

    _record_expected_failure(
        evaluator,
        calls,
        name="explicit_nonzero_gate_false_rejected",
        values=exact_inputs,
        gate=lambda _L: False,
        expected_exception=(ValueError,),
    )
    _record_expected_failure(
        evaluator,
        calls,
        name="noncallable_gate_rejected",
        values=exact_inputs,
        gate=1,
        expected_exception=(TypeError,),
    )
    _record_expected_failure(
        evaluator,
        calls,
        name="nonliteral_gate_result_rejected",
        values=exact_inputs,
        gate=lambda _L: 1,
        expected_exception=(ValueError,),
    )

    zero_denominator_inputs = list(exact_inputs)
    zero_denominator_inputs[0] = Fraction(1)
    zero_denominator_inputs[15] = Fraction(1)  # L = 1 - delta0*Z^2 = 0.
    zero_denominator_inputs[12] = Fraction(1)
    _record_expected_failure(
        evaluator,
        calls,
        name="zero_denominator_requires_false_gate",
        values=tuple(zero_denominator_inputs),
        gate=lambda L: L != 0,
        expected_exception=(ValueError,),
    )

    interval_zero_inputs = list(interval_inputs)
    interval_zero_inputs[12] = RationalInterval.around(Fraction(1), Fraction(1, 1000))
    interval_zero_inputs[15] = RationalInterval.around(Fraction(1), Fraction(1, 1000))
    _record_expected_failure(
        evaluator,
        calls,
        name="interval_denominator_may_contain_zero_requires_false_gate",
        values=tuple(interval_zero_inputs),
        gate=lambda L: L.lower > 0 or L.upper < 0,
        expected_exception=(ValueError,),
    )

    _record_expected_failure(
        evaluator,
        calls,
        name="delta_y_zero_premise_must_be_explicit",
        values=exact_inputs,
        gate=lambda _L: True,
        delta_y_zero_premise=False,
        expected_exception=(ValueError,),
    )

    float_inputs = list(exact_inputs)
    float_inputs[0] = 1.0
    _record_expected_failure(
        evaluator,
        calls,
        name="builtin_float_rejected",
        values=tuple(float_inputs),
        gate=lambda _L: True,
        expected_exception=(TypeError,),
    )

    exact_tuple = _result_tuple(exact)
    _require(exact_tuple[0] != 0 and exact_tuple[1] != 0, "fixture accidentally zeroed U rows")
    _require(exact_tuple[5] != -exact_tuple[2], "fixture accidentally imposed Q_y = -Q")
    checks = {
        "exact_fraction_values_match_independent_hand_calculation": exact_tuple == HAND_EXPECTED,
        "exact_output_has_seven_fields": len(exact_tuple) == 7,
        "ordinary_m_rows_are_nonzero_and_ordered": exact_inputs[4:7] == (
            Fraction(1), Fraction(3), Fraction(4)
        ),
        "delta_Z_nonzero": exact_inputs[13] != 0,
        "delta_ZZ_nonzero": exact_inputs[14] != 0,
        "gap_U_not_zeroed": exact_inputs[0] != 0 and exact_inputs[1] != 0,
        "gap_Uy_and_my_not_zeroed": all(value != 0 for value in exact_inputs[7:12]),
        "Q_y_not_auto_negated": exact_tuple[5] != -exact_tuple[2],
        "rational_interval_contains_exact_point": calls[1]["passed"],
        "nonzero_gate_false_rejected": calls[2]["passed"],
        "noncallable_gate_rejected": calls[3]["passed"],
        "nonliteral_gate_result_rejected": calls[4]["passed"],
        "zero_denominator_rejected": calls[5]["passed"],
        "interval_zero_denominator_rejected": calls[6]["passed"],
        "delta_y_zero_premise_required": calls[7]["passed"],
        "builtin_float_rejected": calls[8]["passed"],
        "no_builtin_float_in_recorded_fixture": not _contains_builtin_float(calls),
        "no_interval_midpoint_in_recorded_fixture": not _contains_midpoint_field(calls),
    }
    _require(all(checks.values()), "Q-jet fixture contains a failed check")
    return {
        "schema": COMPONENT + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_checks_passed",
        "passed": True,
        "checks_passed": True,
        "synthetic_only": True,
        "coordinate_convention": {
            "coordinate": "xi",
            "held_fixed": "xi",
            "derivative_variable": "Z",
            "derivative_kind": "ordinary fixed-xi Z derivatives",
            "rows_are_divided_taylor_coefficients": False,
        },
        "input_order": {
            "U": ["U0", "U1", "U2"],
            "m": ["m0", "m_Z", "m_ZZ", "m_ZZZ"],
            "U_y": ["Uy0", "Uy1"],
            "m_y": ["my0", "my1", "my2"],
            "delta": ["delta", "delta_Z", "delta_ZZ"],
        },
        "return_order": ["U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ"],
        "gap_zero_assumptions": False,
        "caller_premises": {
            "delta_y_zero_premise_supplied": True,
            "delta_y_zero_proved": False,
            "premise_is_not_source_identity": True,
        },
        "hand_calculation": HAND_CALCULATION,
        "expected_jets": [str(value) for value in HAND_EXPECTED],
        "calls": calls,
        "call_count": len(calls),
        "checks": checks,
        "interval_backend": {
            "type": "independent RationalInterval",
            "exact_rational_endpoints": True,
            "directed_rounding_certified": False,
            "midpoint_used": False,
        },
        "source_identity_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_value_authentication": False,
        "authentication_proved": False,
        "source_graph_traversed": False,
        "ancestor_replayed": False,
        "numeric_use_admissible": False,
        "certified": False,
    }


def _validate_fixture(value: Mapping[str, Any]) -> None:
    _require(value.get("passed") is True and value.get("checks_passed") is True, "fixture did not pass")
    _require(value.get("synthetic_only") is True, "fixture is not synthetic-only")
    _require(value.get("call_count") == 9, "fixture call count changed")
    calls = value.get("calls")
    _require(
        isinstance(calls, list)
        and len(calls) == 9
        and all(isinstance(item, Mapping) and item.get("invoked") is True and item.get("passed") is True for item in calls),
        "fixture call log contains a failed call",
    )
    _require(not _contains_builtin_float(value), "fixture exposed a built-in float")
    _require(not _contains_midpoint_field(value), "fixture exposed an interval midpoint")
    for key in (
        "source_identity_proved",
        "source_values_selected",
        "source_values_computed",
        "source_value_authentication",
        "authentication_proved",
        "source_graph_traversed",
        "ancestor_replayed",
        "numeric_use_admissible",
        "certified",
    ):
        if key in value:
            _require(value[key] is False, "fixture made forbidden claim: " + key)
    checks = value.get("checks")
    _require(
        isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()),
        "fixture contains a failed check",
    )


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTERED,
        {
            "schema": COMPONENT + "_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "new_leaf_imported": False,
            "synthetic_fixture_invocations": 0,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_computed": False,
            "source_identity_proved": False,
            "source_value_authentication": False,
            "runtime_model_metadata": None,
        },
    )

    leaf_imported = False
    fixture_invocations = 0
    calls: list[dict[str, Any]] = []
    fixture_result: Any = None
    try:
        module = _import_leaf()
        leaf_imported = True
        fixture_invocations = 1
        fixture_result = _run_fixture(module, calls)
        _validate_fixture(fixture_result)
        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_imported": leaf_imported,
            "synthetic_fixture_invocations": fixture_invocations,
            "consumer_calls": fixture_invocations,
            "leaf_api_invocations": len(calls),
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "source_value_authentication": False,
            "source_identity_proved": False,
            "numeric_use_admissible": False,
            "certified": False,
            "runtime_model_metadata": None,
            "source_manifest": _json_safe(source),
            "result": _json_safe(fixture_result),
            "error": None,
        }
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_imported": leaf_imported,
            "synthetic_fixture_invocations": fixture_invocations,
            "consumer_calls": fixture_invocations,
            "leaf_api_invocations": len(calls),
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "source_value_authentication": False,
            "source_identity_proved": False,
            "numeric_use_admissible": False,
            "certified": False,
            "runtime_model_metadata": None,
            "source_manifest": _json_safe(source),
            "result": _json_safe(fixture_result),
            "fixture_call_log": _json_safe(calls),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }

    result_sha = _write_exclusive(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "entry_sha256": entry_sha,
                "synthetic_fixture_invocations": fixture_invocations,
                "consumer_calls": fixture_invocations,
                "leaf_api_invocations": len(calls),
                "ancestor_math_calls": 0,
                "source_identity_proved": False,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
