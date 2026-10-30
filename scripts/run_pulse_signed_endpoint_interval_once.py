"""Run the frozen synthetic signed-endpoint interval fixture exactly once.

The runner authenticates exactly two raw inputs before the new leaf is loaded:
the leaf and this runner.  It writes an exclusive ``entered.json`` first,
then calls only the new endpoint evaluator with exact ``Fraction`` values and
one rational-endpoint interval backend.  No old producer, integral, graph,
quadrature, or source value is imported or called.

The coordinated direct interface is ``evaluate_signed_endpoint`` (aliases
are accepted by the adapter).  Its direct inputs are ``R_F``,
``R_F_Z_at_xi``, ``Q_y_F``, ``Q_y_F_Z_at_xi``, ``dQ_y``, ``dQ_y_Z_at_xi``,
``d``, ``d_Z_at_xi`` and caller-given ``expm1_d`` (or an equivalent caller
callback).  The repaired moving-endpoint interface additionally accepts
``xi_e_Z``, ``R_F_xi``, ``Q_y_F_xi``, ``dQ_y_xi`` and ``d_xi``.  The fixture checks

``E = R_F H``
``E_Z|xi = R_F_Z H + R_F H_Z``
``total E_Z = E_Z|xi + xi_e_Z*(R_F_xi H + R_F H_xi)``.

The expm1 values are exact synthetic data supplied by the fixture, not an
evaluation of a real exponential.  A caller certification boolean is kept as
arithmetic metadata only and is never promoted to source authentication,
physical error certification, or source identity.
"""

from __future__ import annotations

import dataclasses
from fractions import Fraction
import hashlib
import importlib.util
import inspect
import json
import re
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "pulse_signed_endpoint_interval_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "pulse_signed_endpoint_interval.py"
RUNNER = ROOT / "scripts" / "run_pulse_signed_endpoint_interval_once.py"
LEAF_REL = "src/openai_ns_reconstruction/pulse_signed_endpoint_interval.py"
RUNNER_REL = "scripts/run_pulse_signed_endpoint_interval_once.py"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class RationalInterval:
    """Closed interval with exact rational endpoints, used only by the fixture."""

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
        reciprocal = (Fraction(1, 1) / other.lower, Fraction(1, 1) / other.upper)
        return self * type(self)(min(reciprocal), max(reciprocal))

    def __rtruediv__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) / self

    def contains(self, value: int | Fraction) -> bool:
        value = Fraction(value)
        return self.lower <= value <= self.upper

    def __repr__(self) -> str:
        return f"RationalInterval({self.lower!r}, {self.upper!r})"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _sha(path: Path) -> str:
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
        )
        + "\n"
    ).encode("utf-8")


def _json_safe(value: Any) -> Any:
    if isinstance(value, Fraction):
        return {"type": "Fraction", "numerator": value.numerator, "denominator": value.denominator}
    if isinstance(value, RationalInterval):
        return {"type": "RationalInterval", "lower": _json_safe(value.lower), "upper": _json_safe(value.upper)}
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _json_safe({field.name: getattr(value, field.name) for field in dataclasses.fields(value)})
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return {str(key): _json_safe(item) for key, item in fields.items()}
    if hasattr(value, "lower") and hasattr(value, "upper"):
        return {"lower": _json_safe(value.lower), "upper": _json_safe(value.upper)}
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
    """Authenticate exactly the new leaf and runner raw bytes."""

    files = source.get("files")
    _require(isinstance(files, list), "source files manifest missing")
    _require(len(files) == 2, "source manifest must contain exactly two files")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        relative = _relative_path(item.get("path"), f"source.files[{index}]")
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            f"source.files[{index}] has no valid SHA-256",
        )
        _require(relative not in hashes, "duplicate source manifest path: " + relative)
        _require(relative in {LEAF_REL, RUNNER_REL}, "manifest contains an out-of-scope file: " + relative)
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
    _require(isinstance(path, str) and isinstance(digest, str), "pending source checkpoint pointer missing")
    return path, digest


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    """Read and authenticate the pending brief/checkpoint before leaf import."""

    _require(not ENTERED.exists() and not RESULT.exists(), "already entered; never re-enter")
    _require(DIRECTORY.is_dir(), "frozen component directory is missing")
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
    for key, expected in (("goal_active", True), ("full_goal_complete", False), ("do_not_reenter", True)):
        _require(brief.get(key) is expected, "pending brief guard failed: " + key)
    _require(source.get("component") == COMPONENT, "source checkpoint component mismatch")
    _require(source.get("one_entry_only") is True, "source checkpoint is not one-entry-only")
    _require(source.get("no_source_graph_traversal") is True, "source graph traversal guard missing")
    if "synthetic_only" in source:
        _require(source.get("synthetic_only") is True, "source checkpoint is not synthetic-only")
    checkpoint_path, checkpoint_expected = _brief_checkpoint(brief)
    checkpoint_relative = _relative_path(checkpoint_path, "brief source checkpoint")
    _require((ROOT / checkpoint_relative).resolve() == SOURCE.resolve(), "brief points to another checkpoint")
    _require(SHA256_RE.fullmatch(checkpoint_expected) is not None, "brief checkpoint SHA malformed")
    source_sha = hashlib.sha256(source_raw).hexdigest()
    _require(source_sha == checkpoint_expected.lower(), "source checkpoint SHA mismatch")
    hashes = _manifest_hashes(source)
    _require(LEAF.is_file(), "new endpoint leaf is missing")
    return dict(source), hashes, source_sha, hashlib.sha256(brief_raw).hexdigest()


def _result_fields(value: Any) -> Mapping[str, Any] | None:
    if isinstance(value, Mapping):
        return value
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {field.name: getattr(value, field.name) for field in dataclasses.fields(value)}
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return fields
    if hasattr(value, "lower") and hasattr(value, "upper"):
        return {"lower": value.lower, "upper": value.upper}
    return None


def _contains_builtin_float(value: Any) -> bool:
    if type(value) is float:
        return True
    if isinstance(value, Mapping):
        return any(_contains_builtin_float(item) for item in value.values())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_builtin_float(item) for item in value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(_contains_builtin_float(getattr(value, field.name)) for field in dataclasses.fields(value))
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return _contains_builtin_float(fields)
    if hasattr(value, "lower") and hasattr(value, "upper"):
        return _contains_builtin_float(value.lower) or _contains_builtin_float(value.upper)
    return False


def _contains_midpoint_field(value: Any) -> bool:
    if isinstance(value, Mapping):
        return any("midpoint" in str(key).lower() or _contains_midpoint_field(item) for key, item in value.items())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_midpoint_field(item) for item in value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(
            "midpoint" in field.name.lower() or _contains_midpoint_field(getattr(value, field.name))
            for field in dataclasses.fields(value)
        )
    fields = getattr(value, "__dict__", None)
    return isinstance(fields, dict) and _contains_midpoint_field(fields)


def _as_fraction(value: Any) -> Fraction | None:
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value, 1)
    return None


def _as_interval(value: Any) -> RationalInterval | None:
    if isinstance(value, RationalInterval):
        return value
    fields = _result_fields(value)
    if fields is None:
        return None
    lower = _as_fraction(fields.get("lower"))
    upper = _as_fraction(fields.get("upper"))
    if lower is None or upper is None:
        return None
    return RationalInterval(lower, upper)


def _field(value: Any, names: tuple[str, ...], label: str) -> Any:
    fields = _result_fields(value)
    _require(fields is not None, "leaf result has no fields: " + label)
    for name in names:
        if name in fields:
            return fields[name]
    raise RuntimeError("leaf result lacks field: " + label)


_RESULT_FIELD_ALIASES = {
    "endpoint": ("endpoint", "E", "value"),
    "fixed_coordinate_partial": (
        "fixed_coordinate_partial",
        "fixed_coordinate_Z1",
        "endpoint_fixed_coordinate_partial",
        "partial",
        "partial_Z",
        "E_Z_fixed",
    ),
    "moving_endpoint_total": (
        "moving_endpoint_total",
        "moving_endpoint_Z1",
        "endpoint_moving_total",
        "endpoint_Z",
        "total",
        "total_Z",
    ),
}


def _given_expm1(value: Any, calls: list[Any]) -> Any:
    """Return exact fixture data; this is not an exponential evaluation."""

    calls.append(value)
    if isinstance(value, RationalInterval):
        if value.contains(Fraction(0)):
            return RationalInterval.around(Fraction(0), Fraction(1, 1000))
        if value.contains(Fraction(1, 2)):
            return RationalInterval.around(Fraction(1, 2), Fraction(1, 1000))
        raise ValueError("fixture expm1 data is undefined for this interval")
    exact = _as_fraction(value)
    if exact == 0:
        return Fraction(0)
    if exact == Fraction(1, 2):
        return Fraction(1, 2)
    raise ValueError("fixture expm1 data is undefined for this d")


class _GivenExpm1:
    """Callable adapter used only if the leaf requests a callback."""

    def __init__(self, calls: list[Any]) -> None:
        self.calls = calls

    def __call__(self, value: Any) -> Any:
        return _given_expm1(value, self.calls)


def _case(
    *,
    R_F: Any,
    R_F_Z: Any,
    Q_y_F: Any,
    Q_y_F_Z: Any,
    dQ_y: Any,
    dQ_y_Z: Any,
    d: Any,
    d_Z: Any,
    expm1_d: Any,
    expm1_calls: list[Any],
    interval_rounding_certified: bool = False,
    moving: bool = False,
) -> dict[str, Any]:
    callback = _GivenExpm1(expm1_calls)
    result = {
        "R_F": R_F,
        "R_F_Z": R_F_Z,
        "Q_y_F": Q_y_F,
        "Q_y_F_Z": Q_y_F_Z,
        "dQ_y": dQ_y,
        "dQ_y_Z": dQ_y_Z,
        "d": d,
        "d_Z": d_Z,
        "expm1_d": expm1_d,
        "expm1": expm1_d,
        "expm1_callback": callback,
        "interval_rounding_certified": interval_rounding_certified,
    }
    if moving:
        result.update(
            {
                "moving_endpoint": True,
                "xi_e_Z": Fraction(2),
                "R_F_xi": Fraction(4),
                "Q_y_F_xi": Fraction(3),
                "dQ_y_xi": Fraction(1),
                "d_xi": Fraction(1, 4),
            }
        )
    return result


def _expected(case: Mapping[str, Any], expm1_d: Fraction) -> dict[str, Fraction]:
    R_F = Fraction(case["R_F"])
    R_F_Z = Fraction(case["R_F_Z"])
    qF = Fraction(case["Q_y_F"])
    qF_Z = Fraction(case["Q_y_F_Z"])
    dQ = Fraction(case["dQ_y"])
    dQ_Z = Fraction(case["dQ_y_Z"])
    d_Z = Fraction(case["d_Z"])
    exp_d = Fraction(1) + expm1_d
    H = dQ + expm1_d * (qF + dQ)
    H_Z = dQ_Z + exp_d * d_Z * (qF + dQ) + expm1_d * (qF_Z + dQ_Z)
    partial = R_F_Z * H + R_F * H_Z
    total = partial
    result = {
        "endpoint": R_F * H,
        "fixed_coordinate_partial": partial,
    }
    if "xi_e_Z" in case:
        H_xi = (
            Fraction(case["dQ_y_xi"])
            + exp_d * Fraction(case["d_xi"]) * (qF + dQ)
            + expm1_d * (Fraction(case["Q_y_F_xi"]) + Fraction(case["dQ_y_xi"]))
        )
        E_xi = Fraction(case["R_F_xi"]) * H + R_F * H_xi
        total += Fraction(case["xi_e_Z"]) * E_xi
        result["moving_endpoint_total"] = total
    return result


def _find_evaluator(leaf: Any) -> Callable[..., Any]:
    for name in (
        "evaluate_signed_endpoint",
        "evaluate_signed_endpoint_interval",
        "evaluate_endpoint",
        "signed_endpoint_and_Z1",
        "ordinary_z1_endpoint",
        "compute_signed_endpoint",
        "pulse_signed_endpoint",
    ):
        candidate = getattr(leaf, name, None)
        if callable(candidate):
            return candidate
    raise RuntimeError("new endpoint leaf has no signed endpoint evaluator")


def _compact(name: str) -> str:
    return name.lower().replace("_", "")


def _call_evaluator(function: Callable[..., Any], case: Mapping[str, Any]) -> Any:
    """Bind the direct endpoint API without retrying a call."""

    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("endpoint evaluator signature is unavailable") from exc
    direct: dict[str, set[str]] = {
        "R_F": {"rf", "radiusf"},
        "R_F_Z": {"rfz", "radiusfz", "rfzatxi", "radiusfzatxi"},
        "Q_y_F": {"qyf", "frameqy"},
        "Q_y_F_Z": {"qyfz", "frameqyz", "qyfzatxi", "frameqyzatxi"},
        "dQ_y": {"dqy", "deltaqy", "differencedqy"},
        "dQ_y_Z": {"dqyz", "deltaqyz", "differencedqyz", "dqyzatxi", "deltaqyzatxi"},
        "d": {"d", "logd", "logradiusratio"},
        "d_Z": {"dz", "logradiusratioz", "dzatxi", "logradiusratiozatxi"},
        "expm1_d": {"expm1d", "expm1value", "expm1dvalue"},
        "expm1": {"expm1"},
        "expm1_callback": {"expm1callback", "expm1fn"},
        "interval_rounding_certified": {"intervalroundingcertified"},
        "moving_endpoint": {"movingendpoint", "moving", "requestmoving"},
        "xi_e_Z": {"xiz", "xiez", "endpointz", "movingz"},
        "R_F_xi": {"rfxi", "radiusfxi"},
        "Q_y_F_xi": {"qyf xi".replace(" ", ""), "qyf_xi".replace("_", "")},
        "dQ_y_xi": {"dqyxi", "deltaqyxi"},
        "d_xi": {"dxi", "logradiusratioxi"},
    }
    case_names = {"case", "data", "payload", "inputs", "endpointdata"}
    kwargs: dict[str, Any] = {}
    consumed: set[str] = set()
    for parameter in signature.parameters.values():
        if parameter.kind == inspect.Parameter.VAR_POSITIONAL:
            raise RuntimeError("endpoint evaluator may not use variadic positional inputs")
        if parameter.kind == inspect.Parameter.POSITIONAL_ONLY:
            raise RuntimeError("endpoint evaluator must expose keyword inputs")
        if parameter.kind == inspect.Parameter.VAR_KEYWORD:
            continue
        name = _compact(parameter.name)
        if name in {_compact(item) for item in case_names}:
            kwargs[parameter.name] = case
            consumed.add("case")
            continue
        matched = None
        for canonical, aliases in direct.items():
            if name == _compact(canonical) or name in aliases:
                matched = canonical
                break
        if matched is not None and matched in case:
            kwargs[parameter.name] = case[matched]
            consumed.add(matched)
            continue
        if parameter.default is not inspect.Parameter.empty:
            continue
        raise RuntimeError("unsupported required endpoint evaluator argument: " + parameter.name)
    _require(kwargs or consumed, "endpoint evaluator declares no coordinated inputs")
    return function(**kwargs)


def _scope_guard(value: Any) -> None:
    fields = _result_fields(value)
    if fields is None:
        return
    for name in ("source_identity_proved", "source_values_bound", "authentication_proved", "source_values_selected", "source_values_computed", "source_graph_traversed", "physical_error_certified", "error_certified", "quadrature_performed", "midpoint_used"):
        if name in fields:
            _require(fields[name] is False, "leaf result made forbidden claim: " + name)
    if "certified" in fields and fields["certified"] is True:
        _require(fields.get("source_identity_proved") is False, "arithmetic certification became source identity")
        _require(fields.get("authentication_proved") is False, "arithmetic certification became authentication")


def _invoke(
    evaluator: Callable[..., Any],
    calls: list[dict[str, Any]],
    *,
    name: str,
    case: Mapping[str, Any],
    expected: Mapping[str, Fraction] | None = None,
    expected_exception: tuple[type[BaseException], ...] | None = None,
    interval_expected: Mapping[str, Fraction] | None = None,
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "invoked": False,
        "expected_exception": [item.__name__ for item in expected_exception] if expected_exception else None,
        "input_scope": "synthetic_fraction_or_rational_interval",
    }
    calls.append(record)
    try:
        record["invoked"] = True
        value = _call_evaluator(evaluator, case)
        record["returned"] = True
        _scope_guard(value)
        _require(not _contains_builtin_float(value), "leaf result exposed a built-in float")
        _require(not _contains_midpoint_field(value), "leaf result exposed an interval midpoint")
        record["output"] = _json_safe(value)
        if expected_exception is not None:
            raise RuntimeError(name + " unexpectedly returned")
        for field_name, expected_value in (expected or {}).items():
            observed = _field(value, _RESULT_FIELD_ALIASES.get(field_name, (field_name,)), field_name)
            observed_fraction = _as_fraction(observed)
            _require(observed_fraction == expected_value, name + " " + field_name + " mismatch")
        if expected is not None:
            record["expected"] = {name: str(item) for name, item in expected.items()}
            record["exact_output"] = True
        if interval_expected is not None:
            for field_name, expected_value in interval_expected.items():
                observed = _as_interval(_field(value, _RESULT_FIELD_ALIASES.get(field_name, (field_name,)), field_name))
                _require(observed is not None, name + " did not return rational interval for " + field_name)
                _require(observed.contains(expected_value), name + " interval misses point for " + field_name)
            record["interval_contains_point"] = True
        record["passed"] = True
        return value
    except BaseException as exc:
        record["returned"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        if expected_exception is not None and isinstance(exc, expected_exception):
            record["passed"] = True
            return None
        record["passed"] = False
        raise


def _interval_case(base: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in base.items():
        if key in {"expm1_d", "expm1"}:
            result[key] = RationalInterval.around(Fraction(value), Fraction(1, 1000))
        elif isinstance(value, (int, Fraction)) and not isinstance(value, bool):
            result[key] = RationalInterval.around(Fraction(value), Fraction(1, 1000))
        else:
            result[key] = value
    result["interval_rounding_certified"] = False
    return result


def _run_fixture(leaf: Any) -> dict[str, Any]:
    evaluator = _find_evaluator(leaf)
    calls: list[dict[str, Any]] = []
    callback_calls: list[Any] = []
    base = _case(
        R_F=Fraction(2), R_F_Z=Fraction(3), Q_y_F=Fraction(5), Q_y_F_Z=Fraction(7),
        dQ_y=Fraction(11), dQ_y_Z=Fraction(13), d=Fraction(0), d_Z=Fraction(0),
        expm1_d=Fraction(0), expm1_calls=callback_calls,
    )
    expected_base = _expected(base, Fraction(0))
    fixed = _invoke(evaluator, calls, name="fixed_xi_d_zero", case=base, expected=expected_base)
    _require(expected_base["endpoint"] == Fraction(22), "fixed hand endpoint changed")
    _require(expected_base["fixed_coordinate_partial"] == Fraction(59), "fixed hand partial changed")
    _require(expected_base["fixed_coordinate_partial"] == Fraction(59), "fixed hand total changed")
    _require(
        _as_fraction(_field(fixed, ("endpoint_Z", "fixed_coordinate_partial"), "fixed endpoint_Z"))
        == Fraction(59),
        "fixed endpoint_Z changed",
    )

    dZ_case = dict(base)
    dZ_case["d_Z"] = Fraction(1, 2)
    dZ_case["expm1_calls"] = callback_calls
    expected_dZ = _expected(dZ_case, Fraction(0))
    dZ_result = _invoke(evaluator, calls, name="fixed_xi_nonzero_d_Z", case=dZ_case, expected=expected_dZ)
    _require(expected_dZ["fixed_coordinate_partial"] == Fraction(75), "d_Z hand partial changed")
    _require(
        _as_fraction(_field(dZ_result, ("endpoint_Z", "fixed_coordinate_partial"), "d_Z endpoint_Z"))
        == Fraction(75),
        "d_Z endpoint_Z changed",
    )

    moving = dict(base)
    moving.update(
        {
            "moving_endpoint": True,
            "xi_e_Z": Fraction(2),
            "R_F_xi": Fraction(4),
            "Q_y_F_xi": Fraction(3),
            "dQ_y_xi": Fraction(1),
            "d_xi": Fraction(1, 4),
        }
    )
    expected_moving = _expected(moving, Fraction(0))
    _invoke(evaluator, calls, name="moving_endpoint_complete", case=moving, expected=expected_moving)
    _require(expected_moving["moving_endpoint_total"] == Fraction(167), "moving hand total changed")

    for missing in ("R_F_xi", "Q_y_F_xi", "dQ_y_xi", "d_xi"):
        missing_case = dict(moving)
        del missing_case[missing]
        _invoke(
            evaluator,
            calls,
            name="moving_missing_" + missing,
            case=missing_case,
            expected_exception=(KeyError, TypeError, ValueError),
        )

    zero = dict(base)
    zero.update({"dQ_y": Fraction(0), "dQ_y_Z": Fraction(0)})
    expected_zero = _expected(zero, Fraction(0))
    _invoke(evaluator, calls, name="zero_signed_difference", case=zero, expected=expected_zero)
    _require(expected_zero["endpoint"] == 0 and expected_zero["fixed_coordinate_partial"] == 0, "zero hand result changed")

    lower_calls: list[Any] = []
    lower = _case(
        R_F=Fraction(2), R_F_Z=Fraction(3), Q_y_F=Fraction(5), Q_y_F_Z=Fraction(7),
        dQ_y=Fraction(11), dQ_y_Z=Fraction(13), d=Fraction(0), d_Z=Fraction(0),
        expm1_d=Fraction(0), expm1_calls=lower_calls,
    )
    upper_calls: list[Any] = []
    upper = _case(
        R_F=Fraction(3), R_F_Z=Fraction(1), Q_y_F=Fraction(4), Q_y_F_Z=Fraction(6),
        dQ_y=Fraction(2), dQ_y_Z=Fraction(5), d=Fraction(0), d_Z=Fraction(0),
        expm1_d=Fraction(0), expm1_calls=upper_calls,
    )
    lower_expected = _expected(lower, Fraction(0))
    upper_expected = _expected(upper, Fraction(0))
    lower_result = _invoke(evaluator, calls, name="signed_interval_lower_endpoint", case=lower, expected=lower_expected)
    upper_result = _invoke(evaluator, calls, name="signed_interval_upper_endpoint", case=upper, expected=upper_expected)
    lower_value = _as_fraction(_field(lower_result, ("endpoint", "E"), "lower endpoint"))
    upper_value = _as_fraction(_field(upper_result, ("endpoint", "E"), "upper endpoint"))
    lower_Z = _as_fraction(_field(lower_result, ("endpoint_Z", "moving_endpoint_total"), "lower endpoint_Z"))
    upper_Z = _as_fraction(_field(upper_result, ("endpoint_Z", "moving_endpoint_total"), "upper endpoint_Z"))
    _require(lower_value is not None and upper_value is not None and lower_Z is not None and upper_Z is not None, "signed endpoints were not exact")
    signed_value = -(upper_value - lower_value)
    signed_Z = -(upper_Z - lower_Z)
    _require(signed_value == Fraction(16) and signed_Z == Fraction(42), "endpoint signs were pre-cancelled or reversed")

    known_d_calls: list[Any] = []
    known_d = _case(
        R_F=Fraction(2), R_F_Z=Fraction(3), Q_y_F=Fraction(5), Q_y_F_Z=Fraction(7),
        dQ_y=Fraction(11), dQ_y_Z=Fraction(13), d=Fraction(1, 2), d_Z=Fraction(1, 2),
        expm1_d=Fraction(1, 2), expm1_calls=known_d_calls,
    )
    expected_known_d = _expected(known_d, Fraction(1, 2))
    _invoke(evaluator, calls, name="given_exact_expm1_nonzero_case", case=known_d, expected=expected_known_d)

    interval_calls: list[Any] = []
    interval = _interval_case(base)
    _invoke(
        evaluator,
        calls,
        name="rational_interval_backend",
        case=interval,
        interval_expected=expected_base,
    )

    certified_calls: list[Any] = []
    certified_case = _case(
        R_F=Fraction(2), R_F_Z=Fraction(3), Q_y_F=Fraction(5), Q_y_F_Z=Fraction(7),
        dQ_y=Fraction(11), dQ_y_Z=Fraction(13), d=Fraction(0), d_Z=Fraction(0),
        expm1_d=Fraction(0), expm1_calls=certified_calls,
        interval_rounding_certified=True,
    )
    certified_result = _invoke(evaluator, calls, name="caller_certification_flag_scope", case=certified_case, expected=expected_base)
    certified_fields = _result_fields(certified_result) or {}
    if "certified" in certified_fields:
        _require(certified_fields["certified"] in (False, True), "certified marker is not boolean")
        _require(certified_fields.get("interval_rounding_certified") is True, "caller certification flag was not retained as metadata")
        _require(certified_fields["certified"] is False, "caller certification flag became a result certificate")
        _require(certified_fields.get("source_identity_proved") is False, "caller bool became source identity")
        _require(certified_fields.get("authentication_proved") is False, "caller bool became authentication")

    given_expm1_inputs = [
        base["expm1_d"],
        known_d["expm1_d"],
        interval["expm1_d"],
        certified_case["expm1_d"],
    ]
    _require(any(_as_fraction(item) == 0 for item in given_expm1_inputs), "d=0 exact expm1 data was not exercised")
    _require(any(_as_fraction(item) == Fraction(1, 2) for item in given_expm1_inputs), "known nonzero expm1 data was not exercised")
    _require(all(not _contains_builtin_float(item) for item in callback_calls + known_d_calls + interval_calls + certified_calls), "fixture callback saw a float")

    checks = {
        "fixed_xi_formula": fixed is not None and expected_base["fixed_coordinate_partial"] == Fraction(59),
        "fixed_xi_nonzero_d_Z": expected_dZ["fixed_coordinate_partial"] == Fraction(75),
        "moving_endpoint_formula": expected_moving["moving_endpoint_total"] == Fraction(167),
        "moving_missing_each_xi_partial_rejected": all(calls[index]["passed"] for index in range(3, 7)),
        "zero_signed_difference": calls[7]["passed"],
        "endpoint_sign_not_pre_cancelled": signed_value == Fraction(16) and signed_Z == Fraction(42),
        "direct_signed_difference_nonzero": expected_base["endpoint"] != 0,
        "nonzero_d_Z": expected_dZ["fixed_coordinate_partial"] != expected_base["fixed_coordinate_partial"],
        "RF_Z_retained": Fraction(3) * Fraction(11) != 0,
        "rational_interval_contains_point": calls[11]["passed"],
        "given_expm1_d_zero": any(_as_fraction(item) == 0 for item in given_expm1_inputs),
        "given_expm1_known_case": any(_as_fraction(item) == Fraction(1, 2) for item in given_expm1_inputs),
        "caller_certification_not_source_cert": calls[12]["passed"],
        "no_builtin_float": not _contains_builtin_float(calls),
        "no_interval_midpoint": not _contains_midpoint_field(calls),
    }
    _require(all(checks.values()), "signed endpoint fixture contains a failed check")
    return {
        "schema": COMPONENT + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_checks_passed",
        "passed": True,
        "checks_passed": True,
        "synthetic_only": True,
        "calls": calls,
        "call_count": len(calls),
        "checks": checks,
        "hand_calculation": {
            "d_zero": {"E": "22", "fixed_coordinate_partial": "59", "moving_total": None},
            "d_Z_half": {"fixed_partial": "75"},
            "moving": {"H": "11", "H_xi": "5", "E_xi": "54", "total": "167"},
            "signed_upper_minus_lower_with_minus_sign": {"E": "16", "E_Z": "42"},
        },
        "expm1_callback_scope": {
            "exact_rational_given_data": True,
            "actual_expm1_evaluated": False,
            "callback_inputs": _json_safe(given_expm1_inputs + callback_calls + known_d_calls + interval_calls + certified_calls),
        },
        "source_identity_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_graph_traversed": False,
        "physical_error_certified": False,
        "quadrature_performed": False,
        "certified": False,
    }


def _validate_fixture(value: Any) -> None:
    _require(isinstance(value, Mapping), "fixture did not return an object")
    _require(value.get("passed") is True and value.get("checks_passed") is True, "fixture checks did not pass")
    _require(value.get("synthetic_only") is True, "fixture is not synthetic-only")
    calls = value.get("calls")
    _require(isinstance(calls, list) and len(calls) == 13, "fixture call log incomplete")
    _require(all(isinstance(item, Mapping) and item.get("invoked") is True and item.get("passed") is True for item in calls), "fixture call log contains a failed call")
    _require(not _contains_builtin_float(value), "fixture exposed a built-in float")
    _require(not _contains_midpoint_field(value), "fixture exposed an interval midpoint")
    for key in ("source_identity_proved", "source_values_selected", "source_values_computed", "source_graph_traversed", "physical_error_certified", "quadrature_performed", "certified"):
        if key in value:
            _require(value[key] is False, "fixture made forbidden claim: " + key)
    checks = value.get("checks")
    _require(isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()), "fixture contains a failed check")


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTERED,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "fixture_invocations": 0,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_computed": False,
            "runtime_model_metadata": None,
        },
    )
    fixture_invocations = 0
    leaf_imported = False
    fixture_result: Any = None
    try:
        spec = importlib.util.spec_from_file_location("pulse_signed_endpoint_interval_new_leaf", LEAF)
        _require(spec is not None and spec.loader is not None, "new endpoint leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        leaf_imported = True
        fixture_invocations = 1
        fixture_result = _run_fixture(leaf)
        _validate_fixture(fixture_result)
        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_imported": leaf_imported,
            "fixture_invocations": fixture_invocations,
            "consumer_calls": fixture_invocations,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "physical_error_certified": False,
            "certified": False,
            "runtime_model_metadata": None,
            "result": fixture_result,
            "error": None,
        }
    except BaseException as exc:
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_imported": leaf_imported,
            "fixture_invocations": fixture_invocations,
            "consumer_calls": fixture_invocations,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "physical_error_certified": False,
            "certified": False,
            "runtime_model_metadata": None,
            "fixture_result": _json_safe(fixture_result),
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        }
    result_sha = _write_exclusive(RESULT, output)
    print(json.dumps({"component": COMPONENT, "phase": PHASE, "status": output["status"], "result_sha256": result_sha, "source_checkpoint_sha256": source_sha, "entry_sha256": entry_sha, "fixture_invocations": fixture_invocations, "source_graph_traversed": False, "certified": False}, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
