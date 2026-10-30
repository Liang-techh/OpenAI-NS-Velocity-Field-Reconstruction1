"""Run the fixed-``y`` F-gap log-scaled leaf exactly once.

This runner is deliberately a synthetic, caller-owned fixture.  It accepts a
pending brief and a frozen source manifest that authenticate exactly this
runner and ``f_gap_fixed_y_log_scaled_jets.py`` by raw SHA-256.  The exclusive
``entered.json`` receipt is written before the leaf is imported.  After that
point the only imported project module is the new leaf, and every invocation
is recorded.  A complete success or failure packet is written once to
``result.json``.

The fixture uses independent exact Fractions.  It checks the third-order
ordinary-Z amplitude and Bell rows, the radial-Q quotient through Q2, the
conditioned Qy rows, fixed-y cancellation in ``R*Qy/Rp``, two independent
nonzero denominator gates, the nonzero-``mu_Z`` fixed-y versus fixed-xi
derivative difference, and rejection of a zero denominator or a built-in
float.  No physical source owner, source identity, physical scale
certification, error deletion, or numerical certification is claimed.

The leaf source is statically parsed with :func:`ast.parse` before its dynamic
import.  This file is intentionally not an executable test invocation in the
worker context; the coordinating root may run it once after freezing the
pending brief, source manifest, and entry directory.
"""

from __future__ import annotations

import ast
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
COMPONENT = "f_gap_fixed_y_log_scaled_jets_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/f_gap_fixed_y_log_scaled_jets.py"
RUNNER_REL = "scripts/run_f_gap_fixed_y_log_scaled_jets_once.py"
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

FIXED_Y = Fraction(22, 1)
SECOND_Y = Fraction(17, 1)
FIXED_XI_REFERENCE = Fraction(11, 1)
Z_VALUE = Fraction(1, 2)
DELTA_ROWS = (Fraction(1, 4), Fraction(0), Fraction(0))
NONZERO_MU_ROWS = (Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16))

BASE_ROWS: dict[str, tuple[Fraction, Fraction, Fraction, Fraction]] = {
    "Pstar": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "u": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "C1": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "C2": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "W11": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "W12": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "log_end_scale": (Fraction(23), Fraction(0), Fraction(0), Fraction(0)),
    "mu": (Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0)),
}

EXPECTED_BASE = {
    "amplitude": (Fraction(-2), Fraction(-8), Fraction(-32), Fraction(-128)),
    "phi": (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    "scaled_m": (Fraction(-2), Fraction(-8), Fraction(-32), Fraction(-128)),
    "scaled_Q": (Fraction(36, 5), Fraction(1784, 75), Fraction(66592, 1125)),
    "scaled_Qy": (Fraction(-36, 5), Fraction(-1784, 75)),
    "log_scale": Fraction(1),
    "log_RQ_scale": Fraction(23),
    "scaled_RQ_y": Fraction(-36, 5),
    "scaled_RQ_y_Z": Fraction(-1784, 75),
}

FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "source_identity",
        "source_identity_proved",
        "source_values_authenticated",
        "source_values_selected",
        "source_values_computed",
        "source_capture_proved",
        "physical_scale",
        "physical_scale_proved",
        "physical_scale_certified",
        "physical_source_owner_proved",
        "source_Rp_owner_authenticated",
        "error_terms_deletable",
        "numeric_use_admissible",
        "certification",
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
    """Closed exact-rational interval used by the synthetic interval call."""

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
    """Authenticate exactly the leaf and runner source files."""

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


def _static_parse_leaf() -> dict[str, Any]:
    """Parse the leaf without importing or executing it."""

    _require(LEAF.is_file(), "new fixed-y leaf is missing")
    raw = LEAF.read_text(encoding="utf-8")
    tree = ast.parse(raw, filename=str(LEAF))
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    _require(
        not any(name.startswith("openai_ns_reconstruction") for name in imported),
        "leaf imports an older project module",
    )
    return {
        "ast_parse_passed": True,
        "leaf_source_bytes": len(raw.encode("utf-8")),
        "leaf_imports": sorted(imported),
        "leaf_imported_before_entry": False,
    }


def _import_leaf() -> Any:
    _require(LEAF.is_file(), "new fixed-y leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_fixed_y_log_scaled_jets_new_leaf", LEAF
    )
    _require(spec is not None and spec.loader is not None,
             "new fixed-y leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new leaf callable is missing: " + name)
    return function


def _result_for_record(value: Any) -> Any:
    as_dict = getattr(value, "as_dict", None)
    if callable(as_dict):
        return as_dict()
    return value


def _invoke(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    *args: Any,
    **kwargs: Any,
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "target": "f_gap_fixed_y_log_scaled_jets.build_f_gap_fixed_y_log_scaled_jets",
        "invoked": True,
    }
    calls.append(record)
    try:
        value = function(*args, **kwargs)
    except BaseException as exc:
        record["returned"] = False
        record["passed"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        raise
    record["returned"] = True
    record["passed"] = True
    record["output"] = _json_safe(_result_for_record(value))
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
        "target": "f_gap_fixed_y_log_scaled_jets.build_f_gap_fixed_y_log_scaled_jets",
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
    record["output"] = _json_safe(_result_for_record(value))
    raise RunnerInputError(name + " unexpectedly returned")


def _rows(value: Any, expected_length: int, label: str) -> tuple[Any, ...]:
    raw = getattr(value, "rows", None)
    if raw is None and isinstance(value, Mapping):
        raw = value.get("rows")
    if raw is None:
        raw = value
    result = tuple(raw)
    _require(len(result) == expected_length,
             f"{label} does not have {expected_length} ordinary rows")
    return result


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "leaf result omitted " + name)
        return value[name]
    result = getattr(value, name, None)
    _require(result is not None, "leaf result omitted " + name)
    return result


def _j4_add(*jets: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(jet[order] for jet in jets) for order in range(4))


def _j4_neg(jet: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(-value for value in jet)


def _j4_product(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
        left[3] * right[0]
        + 3 * left[2] * right[1]
        + 3 * left[1] * right[2]
        + left[0] * right[3],
    )


def _j4_reciprocal(jet: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if jet[0] == 0:
        raise ZeroDivisionError("independent mu reciprocal has zero Z0")
    r0 = Fraction(1) / jet[0]
    r1 = -(jet[1] * r0) / jet[0]
    r2 = -(2 * jet[1] * r1 + jet[2] * r0) / jet[0]
    r3 = -(3 * jet[1] * r2 + 3 * jet[2] * r1 + jet[3] * r0) / jet[0]
    return r0, r1, r2, r3


def _j3_add(*jets: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(jet[order] for jet in jets) for order in range(3))


def _j3_neg(jet: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(-value for value in jet)


def _j3_product(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
    )


def _j3_reciprocal(jet: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if jet[0] == 0:
        raise ZeroDivisionError("independent Q denominator has zero Z0")
    r0 = Fraction(1) / jet[0]
    r1 = -(jet[1] * r0) / jet[0]
    r2 = -(2 * jet[1] * r1 + jet[2] * r0) / jet[0]
    return r0, r1, r2


def _independent_calculation(
    rows: Mapping[str, tuple[Fraction, ...]],
    y: Fraction,
    z: Fraction,
    delta: tuple[Fraction, Fraction, Fraction],
) -> dict[str, Any]:
    coefficient = _j4_add(
        _j4_product(rows["C1"], rows["W11"]),
        _j4_product(rows["C2"], rows["W12"]),
    )
    amplitude = _j4_neg(
        _j4_product(
            _j4_product(rows["Pstar"], rows["u"]),
            coefficient,
        )
    )
    phi = _j4_add(
        rows["log_end_scale"],
        _j4_product(
            (Fraction(13, 2), Fraction(0), Fraction(0), Fraction(0)),
            _j4_reciprocal(rows["mu"]),
        ),
        (Fraction(-13), Fraction(0), Fraction(0), Fraction(0)),
        (-y, Fraction(0), Fraction(0), Fraction(0)),
    )
    phi1, phi2, phi3 = phi[1], phi[2], phi[3]
    bell = (
        Fraction(1),
        phi1,
        phi2 + phi1 * phi1,
        phi3 + 3 * phi1 * phi2 + phi1 * phi1 * phi1,
    )
    scaled_m = _j4_product(amplitude, bell)

    z_jet = (z, Fraction(1), Fraction(0))
    one = (Fraction(1), Fraction(0), Fraction(0))
    z_square = _j3_product(z_jet, z_jet)
    one_minus_z_square = tuple(one[i] - z_square[i] for i in range(3))
    one_minus_delta = tuple(one[i] - delta[i] for i in range(3))
    denominator = tuple(one[i] - _j3_product(delta, z_square)[i] for i in range(3))
    m_value = scaled_m[:3]
    m_z = scaled_m[1:]
    numerator = _j3_neg(
        _j3_add(
            _j3_product(_j3_product(one_minus_delta, z_jet), m_value),
            _j3_product(one_minus_z_square, m_z),
        )
    )
    scaled_Q = _j3_product(numerator, _j3_reciprocal(denominator))
    scaled_Qy = (-scaled_Q[0], -scaled_Q[1])
    return {
        "amplitude": amplitude,
        "phi": phi,
        "scaled_m": scaled_m,
        "scaled_Q": scaled_Q,
        "scaled_Qy": scaled_Qy,
        "log_scale": phi[0],
        "log_RQ_scale": phi[0] + y,
        "scaled_RQ_y": scaled_Qy[0],
        "scaled_RQ_y_Z": scaled_Qy[1],
    }


def _independent_fixed_xi_phi_z(
    log_end_scale: tuple[Fraction, ...],
    mu: tuple[Fraction, ...],
    xi: Fraction,
) -> Fraction:
    """Return the fixed-``xi`` phi_Z row for the coordinate comparison."""

    numerator = (Fraction(13, 2) - xi, Fraction(0), Fraction(0), Fraction(0))
    phi = _j4_add(
        log_end_scale,
        _j4_product(numerator, _j4_reciprocal(mu)),
        (Fraction(-13), Fraction(0), Fraction(0), Fraction(0)),
    )
    return phi[1]


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


def _validate_output(
    value: Any,
    expected: Mapping[str, Any],
    *,
    fixed_y: Fraction,
    interval: bool = False,
) -> None:
    _require(_field(value, "coordinate") == "fixed_y", "leaf changed fixed-y convention")
    _require(_field(value, "fixed_y") == fixed_y if not interval else isinstance(_field(value, "fixed_y"), RationalInterval),
             "leaf changed fixed y")
    _require("ordinary" in _field(value, "row_convention") and "Taylor" in _field(value, "row_convention"),
             "leaf changed ordinary-Z row convention")
    _require(_field(value, "pstar_applied_once") is True,
             "leaf did not mark algebraic Pstar use exactly once")
    _require(_field(value, "scale_materialized") is False,
             "leaf materialized the logarithmic scale")
    for name in (
        "source_formula_conditioned",
        "fixed_y_premise",
        "gap_U_zero_conditioned",
        "gap_Uy_zero_conditioned",
        "my_identity_conditioned",
        "delta_y_zero_conditioned",
        "Rp_Z_zero_premise",
        "rq_y_independent_of_fixed_y",
        "two_distinct_y_rqy_equal_conditioned",
        "nonzero_m_Q_rows_retained",
    ):
        _require(_field(value, name) is True, "leaf omitted premise flag: " + name)
    for name in FORBIDDEN_TRUE_KEYS:
        if hasattr(value, name) or (isinstance(value, Mapping) and name in value):
            _require(_field(value, name) is False, "leaf made a forbidden claim: " + name)

    for name in ("amplitude", "phi", "scaled_m"):
        observed = _rows(_field(value, name), 4, name)
        wanted = expected[name]
        if interval:
            _require(all(isinstance(item, RationalInterval) for item in observed),
                     name + " interval rows were narrowed or materialized")
            _require(all(item.contains(wanted[index]) for index, item in enumerate(observed)),
                     name + " interval misses independent exact value")
        else:
            _require(all(type(item) is Fraction for item in observed),
                     name + " is not exact Fraction output")
            _require(observed == wanted, name + " differs from independent calculation")

    for name in ("scaled_Q",):
        observed = _rows(_field(value, name), 3, name)
        wanted = expected[name]
        if interval:
            _require(all(isinstance(item, RationalInterval) for item in observed),
                     name + " interval rows were narrowed or materialized")
            _require(all(item.contains(wanted[index]) for index, item in enumerate(observed)),
                     name + " interval misses independent exact value")
        else:
            _require(observed == wanted, name + " differs from independent calculation")

    observed_qy = tuple(_field(value, "scaled_Qy"))
    if interval:
        _require(all(isinstance(item, RationalInterval) for item in observed_qy),
                 "scaled_Qy interval rows were narrowed or materialized")
        _require(all(item.contains(expected["scaled_Qy"][index])
                     for index, item in enumerate(observed_qy)),
                 "scaled_Qy interval misses independent exact value")
    else:
        _require(observed_qy == expected["scaled_Qy"],
                 "scaled_Qy differs from independent calculation")

    for name in ("log_scale", "log_RQ_scale", "scaled_RQ_y", "scaled_RQ_y_Z"):
        observed = _field(value, name)
        wanted = expected[name]
        if interval:
            _require(isinstance(observed, RationalInterval) and observed.contains(wanted),
                     name + " interval misses independent exact value")
        else:
            _require(observed == wanted, name + " differs from independent calculation")

    _require(not _contains_builtin_float(_result_for_record(value)),
             "leaf returned a built-in float")
    _require(not _contains_midpoint(_result_for_record(value)),
             "leaf exposed an interval midpoint")


def _fixture_rows() -> dict[str, tuple[Fraction, ...]]:
    return {name: tuple(rows) for name, rows in BASE_ROWS.items()}


def _invoke_with_gates(
    calls: list[dict[str, Any]],
    builder: Callable[..., Any],
    name: str,
    rows: Mapping[str, Any],
    y: Any,
    z: Any,
    delta: Any,
) -> Any:
    mu_trace: list[Any] = []
    q_trace: list[Any] = []

    def mu_gate(value: Any) -> bool:
        mu_trace.append(value)
        return _gate_excludes_zero(value)

    def q_gate(value: Any) -> bool:
        q_trace.append(value)
        return _gate_excludes_zero(value)

    result = _invoke(
        calls,
        name,
        builder,
        *(rows[item] for item in INPUT_NAMES),
        y,
        z,
        delta,
        mu_nonzero_gate=mu_gate,
        q_denominator_nonzero_gate=q_gate,
    )
    calls[-1]["mu_gate_trace"] = _json_safe(mu_trace)
    calls[-1]["q_gate_trace"] = _json_safe(q_trace)
    _require(len(mu_trace) == 1 and len(q_trace) == 1,
             name + " did not call each nonzero gate exactly once")
    return result


def _run_fixture(module: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    builder = _callable(module, "build_f_gap_fixed_y_log_scaled_jets")
    rows = _fixture_rows()
    independent = _independent_calculation(rows, FIXED_Y, Z_VALUE, DELTA_ROWS)
    _require(independent == EXPECTED_BASE,
             "independent exact fixed-y calculation differs from pinned hand values")

    base = _invoke_with_gates(
        calls, builder, "fraction_hand_check", rows, FIXED_Y, Z_VALUE, DELTA_ROWS
    )
    _validate_output(base, independent, fixed_y=FIXED_Y)

    second = _invoke_with_gates(
        calls, builder, "different_fixed_y", rows, SECOND_Y, Z_VALUE, DELTA_ROWS
    )
    independent_second = _independent_calculation(rows, SECOND_Y, Z_VALUE, DELTA_ROWS)
    _validate_output(second, independent_second, fixed_y=SECOND_Y)
    _require(_field(base, "scaled_RQ_y") == _field(second, "scaled_RQ_y"),
             "two fixed-y calls changed scaled_RQ_y")
    _require(_field(base, "scaled_RQ_y_Z") == _field(second, "scaled_RQ_y_Z"),
             "two fixed-y calls changed scaled_RQ_y_Z")
    _require(_field(base, "log_RQ_scale") == _field(second, "log_RQ_scale"),
             "two fixed-y calls changed log_RQ_scale")

    # A nonzero mu_Z row distinguishes the fixed-y derivative from the
    # fixed-xi derivative.  The comparison is independently calculated here,
    # not recovered from the leaf output.
    mu_z_rows = dict(rows)
    mu_z_rows["mu"] = NONZERO_MU_ROWS
    independent_mu_z = _independent_calculation(
        mu_z_rows, FIXED_Y, Z_VALUE, DELTA_ROWS
    )
    nonzero_mu_z = _invoke_with_gates(
        calls,
        builder,
        "nonzero_mu_Z_fixed_y_coordinate_check",
        mu_z_rows,
        FIXED_Y,
        Z_VALUE,
        DELTA_ROWS,
    )
    _validate_output(nonzero_mu_z, independent_mu_z, fixed_y=FIXED_Y)
    fixed_y_phi_z = independent_mu_z["phi"][1]
    fixed_xi_phi_z = _independent_fixed_xi_phi_z(
        rows["log_end_scale"], NONZERO_MU_ROWS, FIXED_XI_REFERENCE
    )
    derivative_difference = fixed_y_phi_z - fixed_xi_phi_z
    expected_transport_difference = (
        -FIXED_XI_REFERENCE
        * NONZERO_MU_ROWS[1]
        / (NONZERO_MU_ROWS[0] * NONZERO_MU_ROWS[0])
    )
    calls[-1]["coordinate_derivative_check"] = _json_safe({
        "fixed_y_phi_Z": fixed_y_phi_z,
        "fixed_xi_phi_Z": fixed_xi_phi_z,
        "difference": derivative_difference,
        "expected_difference": expected_transport_difference,
        "xi": FIXED_XI_REFERENCE,
        "mu0": NONZERO_MU_ROWS[0],
        "mu_Z": NONZERO_MU_ROWS[1],
    })
    _require(fixed_y_phi_z == Fraction(-13, 2),
             "fixed-y phi_Z did not equal -13/2")
    _require(fixed_xi_phi_z == Fraction(9, 2),
             "fixed-xi hand phi_Z did not equal 9/2")
    _require(derivative_difference == Fraction(-11),
             "fixed-y/fixed-xi phi_Z difference did not equal -11")
    _require(derivative_difference == expected_transport_difference,
             "coordinate derivative transport identity changed")

    interval_rows = {
        name: tuple(RationalInterval.around(item) for item in rows[name])
        for name in INPUT_NAMES
    }
    interval_delta = tuple(RationalInterval.around(item) for item in DELTA_ROWS)
    interval = _invoke_with_gates(
        calls,
        builder,
        "rational_interval_hand_check",
        interval_rows,
        RationalInterval.around(FIXED_Y),
        RationalInterval.around(Z_VALUE),
        interval_delta,
    )
    _validate_output(interval, independent, fixed_y=FIXED_Y, interval=True)

    zero_mu_rows = dict(rows)
    zero_mu_rows["mu"] = (Fraction(0), Fraction(0), Fraction(0), Fraction(0))
    _expected_failure(
        calls,
        "zero_mu_gate_rejected",
        builder,
        (ValueError,),
        *(zero_mu_rows[item] for item in INPUT_NAMES),
        FIXED_Y,
        Z_VALUE,
        DELTA_ROWS,
        mu_nonzero_gate=_gate_excludes_zero,
        q_denominator_nonzero_gate=_gate_excludes_zero,
    )

    _expected_failure(
        calls,
        "zero_q_denominator_gate_rejected",
        builder,
        (ValueError,),
        *(rows[item] for item in INPUT_NAMES),
        FIXED_Y,
        Fraction(1),
        (Fraction(1), Fraction(0), Fraction(0)),
        mu_nonzero_gate=_gate_excludes_zero,
        q_denominator_nonzero_gate=_gate_excludes_zero,
    )

    float_rows = dict(rows)
    float_rows["Pstar"] = (1.0, Fraction(1), Fraction(1), Fraction(1))
    _expected_failure(
        calls,
        "builtin_float_rejected",
        builder,
        (TypeError,),
        *(float_rows[item] for item in INPUT_NAMES),
        FIXED_Y,
        Z_VALUE,
        DELTA_ROWS,
        mu_nonzero_gate=_gate_excludes_zero,
        q_denominator_nonzero_gate=_gate_excludes_zero,
    )

    checks = {
        "independent_amplitude": tuple(_rows(_field(base, "amplitude"), 4, "amplitude")) == EXPECTED_BASE["amplitude"],
        "independent_scaled_m": tuple(_rows(_field(base, "scaled_m"), 4, "scaled_m")) == EXPECTED_BASE["scaled_m"],
        "independent_scaled_Q": tuple(_rows(_field(base, "scaled_Q"), 3, "scaled_Q")) == EXPECTED_BASE["scaled_Q"],
        "Qy_equals_minus_Q_through_Z1": tuple(_field(base, "scaled_Qy")) == tuple(-item for item in EXPECTED_BASE["scaled_Q"][:2]),
        "two_distinct_y_values_same_RQy": (
            _field(base, "scaled_RQ_y") == _field(second, "scaled_RQ_y")
            and _field(base, "scaled_RQ_y_Z") == _field(second, "scaled_RQ_y_Z")
        ),
        "two_distinct_y_values_same_log_RQ_scale": _field(base, "log_RQ_scale") == _field(second, "log_RQ_scale"),
        "nonzero_mu_Z_fixed_y_phi_Z": fixed_y_phi_z == Fraction(-13, 2),
        "fixed_xi_comparison_phi_Z": fixed_xi_phi_z == Fraction(9, 2),
        "fixed_y_minus_fixed_xi_difference": derivative_difference == Fraction(-11),
        "mu_transport_difference": derivative_difference == expected_transport_difference,
        "rational_interval_contains_exact_point": True,
        "zero_mu_gate_rejected": calls[-3]["passed"],
        "zero_q_gate_rejected": calls[-2]["passed"],
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
        "coordinate": "fixed_y",
        "input_order": list(INPUT_NAMES) + ["fixed_y", "Z", "delta0", "delta1", "delta2"],
        "hand_values": {
            key: [str(item) for item in value] if isinstance(value, tuple) else str(value)
            for key, value in EXPECTED_BASE.items()
        },
        "independent_hand_calculation_used": True,
        "coordinate_derivative_check": {
            "fixed_y_phi_Z": str(fixed_y_phi_z),
            "fixed_xi_phi_Z": str(fixed_xi_phi_z),
            "difference": str(derivative_difference),
            "expected_difference": str(expected_transport_difference),
            "identity": "fixed_y_phi_Z-fixed_xi_phi_Z=-xi*mu_Z/mu0^2",
        },
        "checks": checks,
        "interval_backend": {
            "type": "independent RationalInterval",
            "exact_rational_endpoints": True,
            "midpoint_used": False,
            "no_interval_midpoint": True,
        },
        "premises": {
            "fixed_y": True,
            "xi_equals_mu_times_y_is_not_fixed_xi": True,
            "ordinary_Z_rows_not_divided_Taylor_coefficients": True,
            "gap_U_zero": True,
            "gap_Uy_zero": True,
            "my_equals_minus_m": True,
            "delta_y_zero": True,
            "Rp_Z_zero_for_RQy_identity": True,
            "gates_are_caller_premises_not_source_proofs": True,
        },
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_capture_proved": False,
        "physical_scale_proved": False,
        "physical_source_owner_proved": False,
        "source_Rp_owner_authenticated": False,
        "error_terms_deletable": False,
        "certification": False,
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
    ast_info: Mapping[str, Any],
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
        "physical_scale_certified": False,
        "physical_source_owner_proved": False,
        "source_Rp_owner_authenticated": False,
        "error_terms_deletable": False,
        "source_parameter_independence_proved": False,
        "numeric_use_admissible": False,
        "certification": False,
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
        "ast_parse_before_import": True,
        "leaf_imported": leaf_imported,
        "fixture_invocations": 1 if leaf_imported else 0,
        "adapter_invocations": len(calls),
        "runtime_model_metadata": None,
        "source_manifest": _json_safe(source),
        "ast": _json_safe(ast_info),
        "fixture": _json_safe(fixture),
        "calls": _json_safe(calls),
        "flags": flags,
        **flags,
        "premises": {
            "synthetic_only": True,
            "fixed_y": True,
            "source_owned_F_bound": False,
            "physical_conclusions": False,
            "uncertified": True,
        },
        "error": _json_safe(error),
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    ast_info = _static_parse_leaf()
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
            "ast_parse_passed": True,
            "new_leaf_imported": False,
            "fixture_invocations": 0,
            "adapter_invocations": 0,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "physical_source_owner_proved": False,
            "error_terms_deletable": False,
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
            ast_info=ast_info,
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
            ast_info=ast_info,
            leaf_imported=leaf_imported,
            error={
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        )

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
        "physical_scale_proved": False,
        "certified": False,
    }, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
