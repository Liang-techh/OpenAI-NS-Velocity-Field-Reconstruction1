"""Enter the caller-owned log-scaled original-gap mean fixture once.

Only the new leaf and this runner are frozen inputs.  The runner never reads
an archive or an old result, and imports only the new leaf after the exclusive
entry receipt.  Its fixture supplies synthetic exact ordinary-Z jets, checks
the returned rows against an independent polynomial calculation, and records
premises and uncertified/source-unbound status separately from the algebraic
check.  The large exponential scale is kept symbolic; this runner never
materializes it.

The coordinated leaf API is::

    build_f_gap_native_log_scaled_mean(
        Pstar, u, C1, C2, W11, W12, log_end_scale, mu, xi,
        *, denominator_nonzero_gate,
    )

where every jet is four ordinary fixed-``xi`` Z rows.  The leaf returns
``FGapNativeLogScaledMean`` with ``amplitude``, ``phi``, ``scaled_m``, and
conditional ``scaled_my`` rows.
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
COMPONENT = "f_gap_native_log_scaled_mean_v3"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/f_gap_native_log_scaled_mean_v3.py"
RUNNER_REL = "scripts/run_f_gap_native_log_scaled_mean_v3_once.py"
LEAF = ROOT / LEAF_REL
RUNNER = ROOT / RUNNER_REL
REQUIRED_INPUTS = frozenset({LEAF_REL, RUNNER_REL})
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

INPUT_NAMES = (
    "Pstar",
    "u",
    "C1",
    "C2",
    "W11",
    "W12",
    "log_end_scale",
    "mu",
)
FIXED_XI = Fraction(11, 1)

# These rows are independent fixture data, not source values.  Every supplied
# source row is nonzero in the base case, including all four rows of mu.
BASE_ROWS: dict[str, tuple[Fraction, Fraction, Fraction, Fraction]] = {
    "Pstar": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "u": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "C1": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "C2": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "W11": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "W12": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "log_end_scale": (Fraction(23), Fraction(-5, 2), Fraction(21, 4), Fraction(23, 8)),
    "mu": (Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)),
}

HAND_EXPECTED = {
    "amplitude": (Fraction(-2), Fraction(-8), Fraction(-32), Fraction(-128)),
    "phi": (Fraction(1), Fraction(2), Fraction(3), Fraction(4)),
    "scaled_m": (Fraction(-2), Fraction(-12), Fraction(-78), Fraction(-548)),
    "scaled_my": (Fraction(2), Fraction(12), Fraction(78)),
    "log_scale": Fraction(1),
}

HAND_DERIVATION = {
    "ordinary_product_rows": {
        "p0": "a0*b0",
        "p1": "a1*b0+a0*b1",
        "p2": "a2*b0+2*a1*b1+a0*b2",
        "p3": "a3*b0+3*a2*b1+3*a1*b2+a0*b3",
    },
    "ones_product": ["1", "2", "4", "8"],
    "mu_reciprocal": ["2", "-1", "1/2", "-1/4"],
    "phi": "log_end_scale + (-9/2)*mu^{-1} - 13 = (1,2,3,4)",
    "exponential_scale_bell_rows": {
        "B2": "phi2+phi1^2 = 7",
        "B3": "phi3+3*phi1*phi2+phi1^3 = 30",
    },
    "scaled_m_rows": {
        "m0": "A0",
        "m1": "A1+A0*phi1",
        "m2": "A2+2*A1*phi1+A0*B2",
        "m3": "A3+3*A2*phi1+3*A1*B2+A0*B3",
    },
    "scale_is_unmaterialized": True,
}

FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "source_identity_proved",
        "source_values_authenticated",
        "source_values_selected",
        "source_values_computed",
        "source_capture_proved",
        "physical_scale_proved",
        "physical_Pstar_owner_bound",
        "numeric_use_admissible",
        "certified",
        "source_graph_traversed",
        "old_source_graph_traversed",
        "ancestor_replayed",
        "old_math_called",
        "source_identified",
        "source_owner_bound",
        "physical_conclusion_proved",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


class RationalInterval:
    """Closed exact-rational interval used only by the synthetic fixture."""

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
    def around(cls, value: Fraction, radius: Fraction = Fraction(1, 1000)) -> "RationalInterval":
        value = Fraction(value)
        radius = Fraction(radius)
        if radius < 0:
            raise ValueError("interval radius must be nonnegative")
        return cls(value - radius, value + radius)

    @staticmethod
    def _coerce(value: Any) -> "RationalInterval":
        if isinstance(value, RationalInterval):
            return value
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise TypeError("interval arithmetic accepts exact rationals only")
        return RationalInterval(Fraction(value), Fraction(value))

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
        reciprocal = (Fraction(1) / other.lower, Fraction(1) / other.upper)
        return self * type(self)(min(reciprocal), max(reciprocal))

    def __rtruediv__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) / self

    def contains(self, value: Fraction) -> bool:
        value = Fraction(value)
        return self.lower <= value <= self.upper

    def __repr__(self) -> str:
        return f"RationalInterval({self.lower!r}, {self.upper!r})"


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
    if isinstance(value, Fraction):
        return {"type": "Fraction", "numerator": value.numerator, "denominator": value.denominator}
    if isinstance(value, RationalInterval):
        return {
            "type": "RationalInterval",
            "lower": _json_safe(value.lower),
            "upper": _json_safe(value.upper),
        }
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return {"type": "builtin_float", "repr": repr(value)}
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return value.as_posix()
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _json_safe(getattr(value, field.name))
            for field in dataclasses.fields(value)
        }
    if callable(value):
        return {"type": "callable", "repr": repr(value)}
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
    path.parent.mkdir(parents=True, exist_ok=True)
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
    _require(not relative.is_absolute() and ".." not in relative.parts,
             label + " must be repository-relative")
    key = relative.as_posix()
    candidate = (ROOT / relative).resolve()
    _require(candidate.is_relative_to(ROOT.resolve()) and candidate.is_file(),
             label + " is missing or escapes repository: " + normalized)
    return key, candidate


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate exactly the new leaf and this runner, with no extras."""

    files = source.get("files")
    _require(isinstance(files, list) and len(files) == 2,
             "source manifest must contain exactly two file entries")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_path(item.get("path"), f"source.files[{index}]")
        _require(key in REQUIRED_INPUTS, "manifest contains an out-of-scope file: " + key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        expected = item.get("sha256")
        _require(isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
                 "invalid SHA-256 for " + key)
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        hashes[key] = actual
    _require(frozenset(hashes) == REQUIRED_INPUTS,
             "manifest must contain exactly leaf and runner raw SHA-256 entries")
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path, expected = checkpoint.get("path"), checkpoint.get("sha256")
    _require(isinstance(raw_path, str) and isinstance(expected, str),
             "pending brief source checkpoint is incomplete")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(),
             "pending brief points to a different source checkpoint")
    actual = _sha256_file(SOURCE)
    _require(SHA256_RE.fullmatch(expected) is not None and actual.lower() == expected.lower(),
             "pending brief source checkpoint SHA-256 mismatch")
    return actual


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(),
             "one-entry component already entered")
    _require(DIRECTORY.is_dir(), "component directory is missing")
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
    source_sha = _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("schema") == COMPONENT + "_checkpoint_v1", "source manifest schema mismatch")
    for key, expected in (
        ("one_entry_only", True),
        ("synthetic_only", True),
        ("no_source_graph_traversal", True),
        ("no_source_identification", True),
    ):
        if key in source:
            _require(source[key] is expected, "source manifest guard changed: " + key)
    for key in FORBIDDEN_TRUE_KEYS:
        if key in source:
            _require(source[key] is False, "source manifest made forbidden claim: " + key)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, hashlib.sha256(brief_raw).hexdigest()


def _import_leaf() -> Any:
    _require(LEAF.is_file(), "new log-scaled mean leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_native_log_scaled_mean_v3_new_leaf", LEAF
    )
    _require(spec is not None and spec.loader is not None,
             "new log-scaled mean leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new leaf callable is missing: " + name)
    return function


def _invoke(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    *args: Any,
    **kwargs: Any,
) -> Any:
    record: dict[str, Any] = {"name": name, "invoked": True}
    calls.append(record)
    try:
        value = function(*args, **kwargs)
    except BaseException as exc:
        record["returned"] = False
        record["passed"] = False
        record["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
        }
        raise
    record["returned"] = True
    record["passed"] = True
    record["output"] = _json_safe(value)
    return value


def _expected_failure(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    expected: tuple[type[BaseException], ...],
    *args: Any,
    **kwargs: Any,
) -> None:
    record: dict[str, Any] = {
        "name": name,
        "invoked": True,
        "expected_exception": [item.__name__ for item in expected],
    }
    calls.append(record)
    try:
        value = function(*args, **kwargs)
    except BaseException as exc:
        record["returned"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["passed"] = isinstance(exc, expected)
        _require(record["passed"], name + " raised an unexpected exception")
        return
    record["returned"] = True
    record["passed"] = False
    record["output"] = _json_safe(value)
    raise RunnerInputError(name + " unexpectedly returned")


def _rows(value: Any, label: str) -> tuple[Any, ...]:
    if isinstance(value, Mapping):
        raw = value.get("rows")
    else:
        raw = getattr(value, "rows", None)
    _require(raw is not None, label + " has no rows")
    result = tuple(raw)
    _require(len(result) == 4, label + " does not have four ordinary rows")
    return result


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "leaf result omitted " + name)
        return value[name]
    result = getattr(value, name, None)
    _require(result is not None, "leaf result omitted " + name)
    return result


def _jet_add(*jets: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(jet[order] for jet in jets) for order in range(4))


def _jet_product(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    """Independent ordinary-derivative product, with explicit 2 and 3 rows."""

    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
        left[3] * right[0]
        + 3 * left[2] * right[1]
        + 3 * left[1] * right[2]
        + left[0] * right[3],
    )


def _jet_reciprocal(jet: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if jet[0] == 0:
        raise ZeroDivisionError("independent reciprocal has zero Z0")
    r0 = Fraction(1) / jet[0]
    r1 = -(jet[1] * r0) / jet[0]
    r2 = -(2 * jet[1] * r1 + jet[2] * r0) / jet[0]
    r3 = -(3 * jet[1] * r2 + 3 * jet[2] * r1 + jet[3] * r0) / jet[0]
    return r0, r1, r2, r3


def _jet_scale(jet: tuple[Fraction, ...], scalar: Fraction) -> tuple[Fraction, ...]:
    return tuple(scalar * value for value in jet)


def _independent_calculation(
    rows: Mapping[str, tuple[Fraction, ...]],
    xi: Fraction,
) -> dict[str, tuple[Fraction, ...] | Fraction]:
    coefficient = _jet_add(
        _jet_product(rows["C1"], rows["W11"]),
        _jet_product(rows["C2"], rows["W12"]),
    )
    amplitude = _jet_scale(
        _jet_product(_jet_product(rows["Pstar"], rows["u"]), coefficient),
        Fraction(-1),
    )
    numerator = (Fraction(13, 2) - xi, Fraction(0), Fraction(0), Fraction(0))
    phi = _jet_add(
        rows["log_end_scale"],
        _jet_product(numerator, _jet_reciprocal(rows["mu"])),
        (Fraction(-13), Fraction(0), Fraction(0), Fraction(0)),
    )
    phi1, phi2, phi3 = phi[1], phi[2], phi[3]
    bell2 = phi2 + phi1 * phi1
    bell3 = phi3 + 3 * phi1 * phi2 + phi1 * phi1 * phi1
    scaled_m = (
        amplitude[0],
        amplitude[1] + amplitude[0] * phi1,
        amplitude[2] + 2 * amplitude[1] * phi1 + amplitude[0] * bell2,
        amplitude[3]
        + 3 * amplitude[2] * phi1
        + 3 * amplitude[1] * bell2
        + amplitude[0] * bell3,
    )
    return {
        "amplitude": amplitude,
        "phi": phi,
        "scaled_m": scaled_m,
        "scaled_my": tuple(-value for value in scaled_m[:3]),
        "log_scale": phi[0],
    }


def _gate_excludes_zero(value: Any) -> bool:
    if isinstance(value, RationalInterval):
        return bool(value.lower > 0 or value.upper < 0)
    if isinstance(value, Fraction):
        return bool(value != 0)
    return False


def _contains_builtin_float(value: Any) -> bool:
    if type(value) is float:
        return True
    if isinstance(value, RationalInterval):
        return False
    if isinstance(value, Mapping):
        return any(_contains_builtin_float(item) for item in value.values())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_builtin_float(item) for item in value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(_contains_builtin_float(getattr(value, field.name))
                   for field in dataclasses.fields(value))
    return False


def _contains_midpoint(value: Any) -> bool:
    if isinstance(value, Mapping):
        for key, item in value.items():
            text = str(key).lower()
            if text in {"no_interval_midpoint", "interval_endpoints_are_exact_and_no_midpoint"}:
                if item is not True:
                    return True
                continue
            if text == "midpoint_used":
                if item is not False:
                    return True
                continue
            if "midpoint" in text and text not in {"no_interval_midpoint"}:
                return True
            if _contains_midpoint(item):
                return True
        return False
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_midpoint(item) for item in value)
    if isinstance(value, RationalInterval):
        return False
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _contains_midpoint({field.name: getattr(value, field.name)
                                   for field in dataclasses.fields(value)})
    return False


def _walk_forbidden_true(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                found.append(path + "." + key_text)
            found.extend(_walk_forbidden_true(item, path + "." + key_text))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_forbidden_true(item, path + "[" + str(index) + "]"))
    return found


def _input_rows_from_hand(hand: Mapping[str, Any]) -> dict[str, tuple[Fraction, ...]]:
    raw_inputs = hand.get("inputs")
    _require(isinstance(raw_inputs, Mapping), "hand-check fixture omitted inputs")
    result: dict[str, tuple[Fraction, ...]] = {}
    for name in INPUT_NAMES:
        rows = _rows(raw_inputs.get(name), "hand-check " + name)
        _require(all(type(value) is Fraction for value in rows),
                 "hand-check " + name + " is not exact rational data")
        _require(tuple(rows) == BASE_ROWS[name],
                 "hand-check " + name + " changed independently pinned fixture data")
        result[name] = tuple(rows)
    _require(raw_inputs.get("xi") == FIXED_XI,
             "hand-check xi is not the fixed synthetic coordinate")
    return result


def _validate_result_packet(
    value: Any,
    expected: Mapping[str, Any],
    *,
    interval: bool = False,
) -> None:
    _require(_field(value, "coordinate") == "fixed_xi", "leaf changed fixed-xi convention")
    _require(_field(value, "pstar_applied_once") is True,
             "leaf did not mark algebraic Pstar use exactly once")
    _require(_field(value, "scale_materialized") is False,
             "leaf materialized the logarithmic scale")
    statement = _field(value, "actual_m_statement")
    _require(
        isinstance(statement, str)
        and "log_scale" in statement
        and "scaled_m[k]" in statement
        and "not materialized" in statement,
        "leaf changed the unmaterialized m statement",
    )
    for name in (
        "source_formula_conditioned",
        "gap_U_zero_conditioned",
        "gap_Uy_zero_conditioned",
        "my_identity_conditioned",
    ):
        _require(_field(value, name) is True, "leaf omitted premise flag: " + name)
    for name in (
        "source_identity_proved",
        "source_values_authenticated",
    ):
        _require(_field(value, name) is False, "leaf made a source claim: " + name)
    _require(_field(value, "xi") == FIXED_XI, "leaf changed fixed xi")
    for name in ("amplitude", "phi", "scaled_m"):
        observed = _rows(_field(value, name), name)
        wanted = expected[name]
        if interval:
            _require(all(isinstance(item, RationalInterval) for item in observed),
                     name + " interval rows were narrowed or materialized")
            _require(all(item.contains(wanted[index]) for index, item in enumerate(observed)),
                     name + " interval misses the independent exact value")
        else:
            _require(all(type(item) is Fraction for item in observed),
                     name + " is not exact Fraction output")
            _require(observed == wanted, name + " differs from independent calculation")
    observed_my = tuple(_field(value, "scaled_my"))
    wanted_my = expected["scaled_my"]
    if interval:
        _require(all(isinstance(item, RationalInterval) for item in observed_my),
                 "scaled_my interval rows were narrowed or materialized")
        _require(all(item.contains(wanted_my[index]) for index, item in enumerate(observed_my)),
                 "scaled_my interval misses the independent exact value")
        log_scale = _field(value, "log_scale")
        _require(isinstance(log_scale, RationalInterval) and log_scale.contains(expected["log_scale"]),
                 "interval log_scale misses phi0")
    else:
        _require(observed_my == wanted_my, "scaled_my differs from independent calculation")
        _require(_field(value, "log_scale") == expected["log_scale"],
                 "log_scale differs from phi0")
    _require(not _contains_builtin_float(value), "leaf returned a built-in float")
    _require(not _contains_midpoint(value), "leaf exposed an interval midpoint")


def _run_fixture(module: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    hand_check = _callable(module, "hand_check_fixture")
    hand = _invoke(calls, "hand_check_fixture", hand_check)
    _require(isinstance(hand, Mapping), "hand-check fixture did not return a mapping")
    rows = _input_rows_from_hand(hand)
    independent = _independent_calculation(rows, FIXED_XI)
    _require(independent == HAND_EXPECTED,
             "independent hand calculation does not match its pinned constants")
    _require(_jet_product(BASE_ROWS["u"], BASE_ROWS["C1"]) ==
             (Fraction(1), Fraction(2), Fraction(4), Fraction(8)),
             "ordinary third-order product coefficients changed")
    _require(all(value != 0 for value in rows["mu"]), "base mu jet contains zero")

    builder = _callable(module, "build_f_gap_native_log_scaled_mean")
    base_args = tuple(rows[name] for name in INPUT_NAMES)
    base = _invoke(
        calls,
        "build_exact_base",
        builder,
        *base_args,
        FIXED_XI,
        denominator_nonzero_gate=_gate_excludes_zero,
    )
    _validate_result_packet(base, independent)

    # Scale every Pstar jet row together.  This isolates one multiplicative
    # Pstar application even though the synthetic base Pstar has nonzero Z
    # rows, rather than changing the Pstar jet shape in the comparison.
    pstar_three = (Fraction(3), Fraction(3), Fraction(3), Fraction(3))
    scaled_args = (pstar_three,) + base_args[1:]
    scaled_independent = _independent_calculation(
        {**rows, "Pstar": pstar_three}, FIXED_XI
    )
    scaled = _invoke(
        calls,
        "build_pstar_three_comparison",
        builder,
        *scaled_args,
        FIXED_XI,
        denominator_nonzero_gate=_gate_excludes_zero,
    )
    _validate_result_packet(scaled, scaled_independent)
    _require(_rows(_field(scaled, "phi"), "scaled phi") ==
             _rows(_field(base, "phi"), "base phi"),
             "Pstar comparison changed phi")
    for name in ("amplitude", "scaled_m"):
        lhs = _rows(_field(scaled, name), name)
        rhs = _rows(_field(base, name), name)
        _require(lhs == tuple(Fraction(3) * item for item in rhs),
                 "Pstar was not applied exactly once to " + name)
    _require(tuple(_field(scaled, "scaled_my")) ==
             tuple(Fraction(3) * item for item in _field(base, "scaled_my")),
             "Pstar was not applied exactly once to scaled_my")

    interval_args = tuple(
        tuple(RationalInterval.around(value) for value in rows[name])
        for name in INPUT_NAMES
    )
    interval = _invoke(
        calls,
        "build_rational_endpoint_interval",
        builder,
        *interval_args,
        FIXED_XI,
        denominator_nonzero_gate=_gate_excludes_zero,
    )
    _validate_result_packet(interval, independent, interval=True)

    zero_mu_args = list(base_args)
    zero_mu_args[7] = (Fraction(0), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16))
    _expected_failure(
        calls,
        "zero_mu_gate_rejected",
        builder,
        (ValueError,),
        *zero_mu_args,
        FIXED_XI,
        denominator_nonzero_gate=_gate_excludes_zero,
    )

    interval_zero_args = list(interval_args)
    interval_zero_args[7] = (
        RationalInterval(Fraction(-1, 1000), Fraction(1, 1000)),
        *interval_args[7][1:],
    )
    _expected_failure(
        calls,
        "interval_containing_zero_gate_rejected",
        builder,
        (ValueError,),
        *interval_zero_args,
        FIXED_XI,
        denominator_nonzero_gate=_gate_excludes_zero,
    )

    float_args = list(base_args)
    float_args[1] = (1.0, Fraction(1), Fraction(1), Fraction(1))
    _expected_failure(
        calls,
        "builtin_float_rejected",
        builder,
        (TypeError,),
        *float_args,
        FIXED_XI,
        denominator_nonzero_gate=_gate_excludes_zero,
    )

    checks = {
        "independent_phi_rows_nonzero_mu": tuple(independent["phi"]) == HAND_EXPECTED["phi"],
        "independent_third_order_product_coefficients": True,
        "exact_amplitude_rows": tuple(_rows(_field(base, "amplitude"), "base amplitude")) == HAND_EXPECTED["amplitude"],
        "exact_scaled_m_rows": tuple(_rows(_field(base, "scaled_m"), "base scaled_m")) == HAND_EXPECTED["scaled_m"],
        "pstar_scaling_once": True,
        "xi_fixed": _field(base, "coordinate") == "fixed_xi" and _field(base, "xi") == FIXED_XI,
        "scale_unmaterialized": _field(base, "scale_materialized") is False,
        "rational_interval_contains_exact_point": True,
        "interval_endpoints_are_exact_and_no_midpoint": not _contains_midpoint(interval),
        "zero_gate_rejected": calls[-3]["passed"],
        "interval_zero_gate_rejected": calls[-2]["passed"],
        "builtin_float_rejected": calls[-1]["passed"],
    }
    _require(all(checks.values()), "synthetic fixture contains a failed check")
    return {
        "schema": COMPONENT + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_checks_passed",
        "passed": True,
        "checks_passed": True,
        "synthetic_only": True,
        "call_count": len(calls),
        "calls": calls,
        "coordinate": "fixed_xi",
        "input_order": list(INPUT_NAMES) + ["xi"],
        "hand_calculation": HAND_DERIVATION,
        "expected_jets": {key: [str(item) for item in value] if isinstance(value, tuple) else str(value)
                          for key, value in HAND_EXPECTED.items()},
        "leaf_hand_expected_ignored": True,
        "independent_hand_calculation_used": True,
        "checks": checks,
        "interval_backend": {
            "type": "independent RationalInterval",
            "exact_rational_endpoints": True,
            "midpoint_used": False,
            "no_interval_midpoint": True,
        },
        "premises": {
            "xi_fixed": True,
            "ordinary_Z_rows_not_divided_Taylor_coefficients": True,
            "gate_is_caller_premise_not_source_proof": True,
            "gap_formula_conditioned": True,
        },
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_capture_proved": False,
        "physical_scale_proved": False,
        "physical_Pstar_owner_bound": False,
        "source_graph_traversed": False,
        "ancestor_replayed": False,
        "numeric_use_admissible": False,
        "certified": False,
        "uncertified": True,
    }


def _validate_fixture(value: Mapping[str, Any]) -> None:
    _require(value.get("passed") is True and value.get("checks_passed") is True,
             "synthetic fixture did not pass")
    _require(value.get("synthetic_only") is True, "fixture is not synthetic-only")
    checks = value.get("checks")
    _require(isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()),
             "fixture contains a failed check")
    _require(value.get("call_count") == 7, "fixture call count changed")
    _require(not _walk_forbidden_true(value), "fixture made a forbidden source/physical claim")
    _require(not _contains_builtin_float(value), "fixture exposed a built-in float")
    _require(not _contains_midpoint(value), "fixture exposed an interval midpoint")


def _base_output(
    *,
    status: str,
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    entry_sha: str,
    calls: list[dict[str, Any]],
    fixture: Any,
    leaf_imported: bool,
    error: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    flags = {
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_capture_proved": False,
        "physical_scale_proved": False,
        "physical_Pstar_owner_bound": False,
        "source_graph_traversed": False,
        "ancestor_replayed": False,
        "numeric_use_admissible": False,
        "certified": False,
        "old_math_called": False,
        "source_identified": False,
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
        "fixture_invocations": 1 if leaf_imported else 0,
        "adapter_invocations": len(calls),
        "runtime_model_metadata": None,
        "source_manifest": _json_safe(source),
        "fixture": _json_safe(fixture),
        "calls": _json_safe(calls),
        "flags": flags,
        **flags,
        "premises": {
            "synthetic_only": True,
            "xi_fixed": True,
            "source_owned_F_bound": False,
            "physical_conclusions": False,
            "uncertified": True,
        },
        "error": _json_safe(error),
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
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
            "new_leaf_imported": False,
            "fixture_invocations": 0,
            "adapter_invocations": 0,
            "source_identity_proved": False,
            "source_graph_traversed": False,
            "numeric_use_admissible": False,
            "certified": False,
            "runtime_model_metadata": None,
        },
    )

    calls: list[dict[str, Any]] = []
    fixture: Any = None
    leaf_imported = False
    try:
        module = _import_leaf()
        leaf_imported = True
        fixture = _run_fixture(module, calls)
        _validate_fixture(fixture)
        output = _base_output(
            status="returned",
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            brief_sha=brief_sha,
            entry_sha=entry_sha,
            calls=calls,
            fixture=fixture,
            leaf_imported=leaf_imported,
        )
    except BaseException as exc:
        output = _base_output(
            status="failed",
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            brief_sha=brief_sha,
            entry_sha=entry_sha,
            calls=calls,
            fixture=fixture,
            leaf_imported=leaf_imported,
            error={
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        )

    # Persist the complete success or exception packet before printing status.
    result_sha = _write_once(RESULT, output)
    print(json.dumps({
        "component": COMPONENT,
        "phase": PHASE,
        "status": output["status"],
        "result_sha256": result_sha,
        "entry_sha256": entry_sha,
        "fixture_invocations": output["fixture_invocations"],
        "adapter_invocations": output["adapter_invocations"],
        "source_identity_proved": False,
        "certified": False,
    }, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
