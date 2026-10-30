"""Run the reduced-history fixed-``y`` F-gap leaf exactly once.

This runner owns a synthetic exact-rational packet only.  The source manifest
must pin exactly this runner and ``f_gap_reduced_history_fixed_y.py`` by raw
SHA-256.  The leaf is parsed with :mod:`ast` before the exclusive entry
receipt is written, and is imported only after that receipt exists.  No
producer, source graph, ancestor, historical error, finite-``N`` object, or
physical F value is read or called.

The main fixture uses ``mu=1/4`` and the non-degenerate synthetic gap
``xi in [11, 12]``.  Since ``y=xi/mu``, its two fixed-``y`` endpoint packets
use ``y=44`` and ``y=48``.  The two endpoint calls share the exact same
caller-owned rows, so equality of their ``R*Q_y/Rp`` rows is only a
conditional same-source algebra check.  The gap integral keeps its decay
factor as an opaque caller-supplied symbol; no ``exp`` or ``expm1`` value is
materialized.

The file is intentionally a one-entry runner, not a worker-side invocation.
The coordinating root freezes the brief, source manifest, and entry directory
before running it once.
"""

from __future__ import annotations

import ast
import dataclasses
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
import re
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_reduced_history_fixed_y_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/f_gap_reduced_history_fixed_y.py"
RUNNER_REL = "scripts/run_f_gap_reduced_history_fixed_y_once.py"
LEAF = ROOT / LEAF_REL
RUNNER = ROOT / RUNNER_REL
REQUIRED_INPUTS = frozenset({LEAF_REL, RUNNER_REL})
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

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
        "source_parameter_independence_proved",
        "source_mu_independence_proved",
        "source_delta_independence_proved",
        "source_Rp_independence_proved",
        "error_certified",
        "error_terms_deletable",
        "numeric_use_admissible",
        "certification",
        "certified",
        "ancestor_replayed",
        "source_graph_traversed",
        "old_source_graph_traversed",
        "old_math_called",
        "physical_conclusion_proved",
        "finite_N_used",
        "finite_N_materialized",
        "historical_error_removed",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


class SymbolicScalar:
    """A tiny exact polynomial in one caller-supplied symbol.

    The signed gap functional accepts a directed ``gap_decay`` input.  The
    fixture passes ``g=1-exp(-4)`` only as a symbolic label; it must never be
    replaced by an arbitrary rational or a materialized exponential.  The
    polynomial representation is enough to independently compare the leaf's
    ordinary-Z1 formulas.
    """

    __slots__ = ("_terms", "label")

    def __init__(self, value: Fraction | int = Fraction(0), *, label: str | None = None) -> None:
        if label is not None:
            self._terms = ((1, Fraction(1)),)
            self.label = label
        else:
            self._terms = ((0, Fraction(value)),) if Fraction(value) else ()
            self.label = None

    @classmethod
    def symbol(cls, label: str) -> "SymbolicScalar":
        return cls(label=label)

    @classmethod
    def _from_terms(cls, terms: Mapping[int, Fraction], label: str | None) -> "SymbolicScalar":
        obj = cls.__new__(cls)
        obj._terms = tuple(sorted((power, Fraction(value)) for power, value in terms.items() if value))
        obj.label = label
        return obj

    @staticmethod
    def _coerce(value: Any) -> "SymbolicScalar":
        if isinstance(value, SymbolicScalar):
            return value
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise TypeError("symbolic arithmetic requires exact rational values")
        return SymbolicScalar(Fraction(value))

    def _combine(self, other: Any, sign: int) -> "SymbolicScalar":
        right = self._coerce(other)
        terms: dict[int, Fraction] = {}
        for power, value in self._terms:
            terms[power] = terms.get(power, Fraction(0)) + value
        for power, value in right._terms:
            terms[power] = terms.get(power, Fraction(0)) + sign * value
        label = self.label or right.label
        return self._from_terms(terms, label)

    def __add__(self, other: Any) -> "SymbolicScalar":
        return self._combine(other, 1)

    def __radd__(self, other: Any) -> "SymbolicScalar":
        return self._coerce(other)._combine(self, 1)

    def __sub__(self, other: Any) -> "SymbolicScalar":
        return self._combine(other, -1)

    def __rsub__(self, other: Any) -> "SymbolicScalar":
        return self._coerce(other)._combine(self, -1)

    def __neg__(self) -> "SymbolicScalar":
        return self._from_terms({power: -value for power, value in self._terms}, self.label)

    def __mul__(self, other: Any) -> "SymbolicScalar":
        right = self._coerce(other)
        terms: dict[int, Fraction] = {}
        for left_power, left_value in self._terms:
            for right_power, right_value in right._terms:
                power = left_power + right_power
                terms[power] = terms.get(power, Fraction(0)) + left_value * right_value
        return self._from_terms(terms, self.label or right.label)

    def __rmul__(self, other: Any) -> "SymbolicScalar":
        return self * other

    def __truediv__(self, other: Any) -> "SymbolicScalar":
        right = self._coerce(other)
        if len(right._terms) != 1 or right._terms[0][0] != 0:
            raise ZeroDivisionError("symbolic denominator is not a nonzero rational")
        denominator = right._terms[0][1]
        if denominator == 0:
            raise ZeroDivisionError("symbolic denominator is zero")
        return self._from_terms(
            {power: value / denominator for power, value in self._terms}, self.label
        )

    def __rtruediv__(self, other: Any) -> "SymbolicScalar":
        return self._coerce(other) / self

    def __eq__(self, other: object) -> bool:
        return isinstance(other, SymbolicScalar) and self._terms == other._terms

    def __repr__(self) -> str:
        if not self._terms:
            return "0"
        pieces: list[str] = []
        for power, coefficient in self._terms:
            if power == 0:
                pieces.append(str(coefficient))
            elif power == 1:
                pieces.append(f"({coefficient})*{self.label or 't'}")
            else:
                pieces.append(f"({coefficient})*{self.label or 't'}^{power}")
        return " + ".join(pieces)


Jet4 = tuple[Fraction, Fraction, Fraction, Fraction]
Jet3 = tuple[Fraction, Fraction, Fraction]

ZERO4: Jet4 = (Fraction(0),) * 4
ONE4: Jet4 = (Fraction(1), Fraction(0), Fraction(0), Fraction(0))
ZERO3: Jet3 = (Fraction(0),) * 3
ONE3: Jet3 = (Fraction(1), Fraction(0), Fraction(0))

PSTAR: Jet4 = (Fraction(1),) * 4
U: Jet4 = (Fraction(1),) * 4
MU: Jet4 = (Fraction(1, 4), Fraction(0), Fraction(0), Fraction(0))
INCOMING_FACTOR1: Jet4 = (Fraction(1),) * 4
INCOMING_M1: Jet4 = (Fraction(1),) * 4
AP: Jet4 = (Fraction(1),) * 4
SCALED_MAIN_ROWS1: Jet4 = (Fraction(1),) * 4
LOG_END_SCALE: Jet4 = ZERO4
Z_VALUE = Fraction(1, 2)
DELTA: Jet3 = (Fraction(1, 4), Fraction(0), Fraction(0))
Y_LOWER = Fraction(44)
Y_UPPER = Fraction(48)
XI_LOWER = Fraction(11)
XI_UPPER = Fraction(12)


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
        return {
            "type": "Fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
        }
    if isinstance(value, SymbolicScalar):
        return {"type": "caller_symbolic", "label": value.label, "expr": repr(value)}
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if type(value) is float:
        return {"type": "builtin_float", "repr": repr(value)}
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _json_safe(getattr(value, field.name))
            for field in dataclasses.fields(value)
        }
    fields = getattr(value, "__dict__", None)
    if isinstance(fields, dict):
        return _json_safe(fields)
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


def _write_exclusive(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return hashlib.sha256(raw).hexdigest()


def _write_atomic(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        _require(not path.exists(), "result.json appeared before atomic commit")
        os.replace(temporary, path)
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    return hashlib.sha256(raw).hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_bytes())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _relative_path(raw_path: Any, label: str) -> str:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    candidate = Path(normalized)
    if candidate.is_absolute():
        try:
            candidate = candidate.resolve().relative_to(ROOT.resolve())
        except ValueError as exc:
            raise RunnerInputError(label + " path is outside repository") from exc
    _require(".." not in candidate.parts, label + " path escapes repository")
    return candidate.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list) and len(files) == 2, "source manifest must contain exactly two files")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        relative = _relative_path(item.get("path"), f"source.files[{index}]")
        _require(relative in REQUIRED_INPUTS, "manifest contains an out-of-scope file: " + relative)
        _require(relative not in hashes, "duplicate manifest path: " + relative)
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + relative,
        )
        path = ROOT / relative
        _require(path.is_file(), "manifest input is missing: " + relative)
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + relative)
        hashes[relative] = actual
    _require(frozenset(hashes) == REQUIRED_INPUTS, "manifest must be exactly leaf and runner")
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        path = brief.get("source_checkpoint_path")
        expected = brief.get("source_checkpoint_sha256")
    else:
        path = checkpoint.get("path")
        expected = checkpoint.get("sha256")
    _require(isinstance(path, str) and isinstance(expected, str), "pending source checkpoint is incomplete")
    candidate = Path(path.replace("\\", "/"))
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(), "pending brief points to a different source checkpoint")
    actual = _sha256_file(SOURCE)
    _require(SHA256_RE.fullmatch(expected) is not None and actual.lower() == expected.lower(), "pending source checkpoint SHA-256 mismatch")
    return actual


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str, dict[str, Any]]:
    _require(not ENTERED.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(DIRECTORY.is_dir(), "component directory is missing")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = _load_object(BRIEF, "pending brief")
    source = _load_object(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    if "status" in brief:
        _require(brief.get("status") in {"pending", "active"}, "pending brief status changed")
    active = brief.get("goal_active") is True or brief.get("active") is True
    _require(active, "pending brief is not active")
    _require(brief.get("full_goal_complete") is False, "pending brief goal is complete")
    _require(brief.get("do_not_reenter") is True, "pending brief do_not_reenter guard missing")
    source_sha = _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("schema") == COMPONENT + "_checkpoint_v1", "source manifest schema mismatch")
    for key in ("one_entry_only", "synthetic_only", "no_source_graph_traversal"):
        _require(source.get(key) is True, "source manifest guard missing or changed: " + key)
    if "no_source_identification" in source:
        _require(source.get("no_source_identification") is True, "source identification guard changed")
    _require(
        source.get("no_ancestor_replay") is True or source.get("no_old_ancestor_replay") is True,
        "source manifest ancestor guard missing or changed",
    )
    for key in FORBIDDEN_TRUE_KEYS:
        if key in source:
            _require(source[key] is False, "source manifest made forbidden claim: " + key)
    hashes = _manifest_hashes(source)
    preflight = {
        "brief_path": BRIEF.as_posix(),
        "brief_phase": brief.get("phase"),
        "brief_active": active,
        "brief_goal_active": brief.get("goal_active"),
        "brief_do_not_reenter": brief.get("do_not_reenter"),
        "source_manifest_path": SOURCE.as_posix(),
        "source_manifest_schema": source.get("schema"),
        "source_checkpoint_sha256": source_sha,
        "brief_sha256": hashlib.sha256(brief_raw).hexdigest(),
        "manifest_paths_exact": sorted(hashes),
        "source_hashes": hashes,
    }
    return source, hashes, source_sha, preflight["brief_sha256"], preflight


def _static_parse_leaf() -> dict[str, Any]:
    _require(LEAF.is_file(), "new reduced-history leaf is missing")
    raw = LEAF.read_text(encoding="utf-8")
    tree = ast.parse(raw, filename=str(LEAF))
    functions: list[dict[str, Any]] = []
    classes: list[str] = []
    assignments: list[str] = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            positional = [item.arg for item in node.args.posonlyargs + node.args.args]
            keyword_only = [item.arg for item in node.args.kwonlyargs]
            functions.append(
                {
                    "name": node.name,
                    "positional": positional,
                    "keyword_only": keyword_only,
                    "defaults": len(node.args.defaults),
                    "keyword_defaults": len(node.args.kw_defaults),
                }
            )
        elif isinstance(node, ast.ClassDef):
            classes.append(node.name)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    assignments.append(target.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            assignments.append(node.target.id)
    function_names = {item["name"] for item in functions}
    _require("build_f_gap_reduced_history_fixed_y" in function_names, "reduced-history builder missing in AST")
    _require("analytic_gap_signed_source_functional" in function_names, "signed functional missing in AST")
    return {
        "parse": "ast.parse",
        "leaf_sha256_at_parse": _sha256_file(LEAF),
        "functions": functions,
        "classes": classes,
        "assignments": assignments,
        "required_api": ["build_f_gap_reduced_history_fixed_y", "analytic_gap_signed_source_functional"],
    }


def _import_leaf() -> Any:
    spec = importlib.util.spec_from_file_location("f_gap_reduced_history_fixed_y_new_leaf", LEAF)
    _require(spec is not None and spec.loader is not None, "new reduced-history leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new leaf callable is missing: " + name)
    return function


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "leaf result omitted " + name)
        return value[name]
    sentinel = object()
    result = getattr(value, name, sentinel)
    _require(result is not sentinel, "leaf result omitted " + name)
    return result


def _rows(value: Any, label: str, length: int) -> tuple[Any, ...]:
    raw = value.get("rows") if isinstance(value, Mapping) else getattr(value, "rows", None)
    _require(raw is not None, label + " has no rows")
    result = tuple(raw)
    _require(len(result) == length, label + " has unexpected row count")
    return result


def _invoke(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    *args: Any,
    **kwargs: Any,
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "invoked": True,
        "returned": False,
        "passed": False,
    }
    calls.append(record)
    try:
        value = function(*args, **kwargs)
        record["returned"] = True
        record["output"] = _json_safe(value)
        record["passed"] = True
        return value
    except BaseException as exc:
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["traceback"] = traceback.format_exc()
        raise


def _expected_failure(
    calls: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    expected_exception: tuple[type[BaseException], ...],
    *args: Any,
    **kwargs: Any,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "name": name,
        "invoked": True,
        "returned": False,
        "passed": False,
        "expected_exception": [item.__name__ for item in expected_exception],
    }
    calls.append(record)
    try:
        value = function(*args, **kwargs)
    except BaseException as exc:
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["traceback"] = traceback.format_exc()
        record["passed"] = isinstance(exc, expected_exception)
        if not record["passed"]:
            raise
        return record
    record["returned"] = True
    record["output"] = _json_safe(value)
    raise RunnerInputError(name + " unexpectedly returned")


def _add4(left: Jet4, right: Jet4) -> Jet4:
    return tuple(left[index] + right[index] for index in range(4))  # type: ignore[return-value]


def _neg4(value: Jet4) -> Jet4:
    return tuple(-item for item in value)  # type: ignore[return-value]


def _scale4(value: Jet4, scalar: Fraction) -> Jet4:
    return tuple(scalar * item for item in value)  # type: ignore[return-value]


def _prod4(left: Jet4, right: Jet4) -> Jet4:
    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
        left[3] * right[0]
        + 3 * left[2] * right[1]
        + 3 * left[1] * right[2]
        + left[0] * right[3],
    )


def _add3(left: Jet3, right: Jet3) -> Jet3:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def _neg3(value: Jet3) -> Jet3:
    return tuple(-item for item in value)  # type: ignore[return-value]


def _prod3(left: Jet3, right: Jet3) -> Jet3:
    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
    )


def _recip4(value: Jet4) -> Jet4:
    x0, x1, x2, x3 = value
    return (
        1 / x0,
        -x1 / (x0 * x0),
        2 * x1 * x1 / (x0 * x0 * x0) - x2 / (x0 * x0),
        -6 * x1 * x1 * x1 / (x0**4)
        + 6 * x1 * x2 / (x0**3)
        - x3 / (x0 * x0),
    )


def _recip3(value: Jet3) -> Jet3:
    x0, x1, x2 = value
    return (
        1 / x0,
        -x1 / (x0 * x0),
        2 * x1 * x1 / (x0**3) - x2 / (x0 * x0),
    )


def _expected_core(y: Fraction, *, pstar: Jet4 = PSTAR, incoming_m1: Jet4 = INCOMING_M1) -> dict[str, Any]:
    incoming_product = _prod4(INCOMING_FACTOR1, incoming_m1)
    selected_product = _prod4(AP, SCALED_MAIN_ROWS1)
    reduced_history = _prod4(MU, _add4(incoming_product, selected_product))
    source_row_target = _neg4(reduced_history)
    amplitude = _prod4(_prod4(pstar, U), reduced_history)
    k = _add4(LOG_END_SCALE, _scale4(_recip4(MU), Fraction(13, 2)))
    k = _add4(k, (Fraction(-13), Fraction(0), Fraction(0), Fraction(0)))
    phi = _add4(k, (-y, Fraction(0), Fraction(0), Fraction(0)))
    phi1, phi2, phi3 = phi[1], phi[2], phi[3]
    bell: Jet4 = (
        Fraction(1),
        phi1,
        phi2 + phi1 * phi1,
        phi3 + 3 * phi1 * phi2 + phi1 * phi1 * phi1,
    )
    scaled_m = _prod4(amplitude, bell)
    z = (Z_VALUE, Fraction(1), Fraction(0))
    z_square = _prod3(z, z)
    one_minus_delta = _add3(ONE3, _neg3(DELTA))
    one_minus_z_square = _add3(ONE3, _neg3(z_square))
    denominator = _add3(ONE3, _neg3(_prod3(DELTA, z_square)))
    m_value = scaled_m[:3]
    m_z = scaled_m[1:]
    numerator = _neg3(_add3(_prod3(_prod3(one_minus_delta, z), m_value), _prod3(one_minus_z_square, m_z)))
    scaled_Q = _prod3(numerator, _recip3(denominator))
    scaled_Qy = _neg3(scaled_Q)[:2]
    return {
        "incoming_product": incoming_product,
        "selected_product": selected_product,
        "reduced_history": reduced_history,
        "source_row_target": source_row_target,
        "amplitude": amplitude,
        "K": k,
        "phi": phi,
        "scaled_m": scaled_m,
        "scaled_Q": scaled_Q,
        "scaled_Qy": scaled_Qy,
        "log_scale": phi[0],
        "log_RQ_scale": k[0],
        "scaled_RQ_y": scaled_Qy[0],
        "scaled_RQ_y_Z": scaled_Qy[1],
    }


def _assert_exact_rows(value: Any, expected: tuple[Fraction, ...], label: str) -> None:
    observed = tuple(value)
    _require(len(observed) == len(expected), label + " row count changed")
    _require(all(type(item) is Fraction for item in observed), label + " is not exact Fraction data")
    _require(observed == expected, label + " differs from independent hand calculation")


def _assert_uncertified(value: Any, label: str) -> None:
    sentinel = object()
    for name in (
        "source_identity",
        "source_identity_proved",
        "source_values_authenticated",
        "physical_scale",
        "physical_scale_proved",
        "physical_scale_certified",
        "physical_source_owner_proved",
        "source_Rp_owner_authenticated",
        "certification",
        "certified",
        "error_certified",
        "error_terms_deletable",
    ):
        observed = value.get(name, sentinel) if isinstance(value, Mapping) else getattr(value, name, sentinel)
        if observed is not sentinel:
            _require(observed is False, label + " made forbidden claim: " + name)


def _gate(return_value: Any = True) -> tuple[Callable[[Any], Any], dict[str, Any]]:
    state: dict[str, Any] = {"calls": 0, "values": []}

    def callback(value: Any) -> Any:
        state["calls"] += 1
        state["values"].append(value)
        return return_value(value) if callable(return_value) else return_value

    return callback, state


def _builder_kwargs(y: Fraction, mu: Jet4 = MU, pstar: Jet4 = PSTAR, incoming_m1: Jet4 = INCOMING_M1) -> tuple[tuple[Any, ...], dict[str, Any]]:
    args: tuple[Any, ...] = (
        pstar,
        U,
        mu,
        INCOMING_FACTOR1,
        incoming_m1,
        AP,
        SCALED_MAIN_ROWS1,
        LOG_END_SCALE,
        Z_VALUE,
        DELTA[0],
        DELTA[1],
        DELTA[2],
        y,
    )
    mu_gate, mu_state = _gate(True)
    q_gate, q_state = _gate(True)
    return args, {
        "mu_nonzero_gate": mu_gate,
        "q_denominator_nonzero_gate": q_gate,
        "_mu_state": mu_state,
        "_q_state": q_state,
    }


def _run_fixture(module: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    builder = _callable(module, "build_f_gap_reduced_history_fixed_y")
    signed = _callable(module, "analytic_gap_signed_source_functional")
    expected_lower = _expected_core(Y_LOWER)
    expected_upper = _expected_core(Y_UPPER)
    _require(XI_LOWER == MU[0] * Y_LOWER and XI_UPPER == MU[0] * Y_UPPER, "xi/mu coordinate chain changed")
    _require(XI_LOWER < XI_UPPER, "synthetic gap is degenerate")

    lower_args, lower_kwargs = _builder_kwargs(Y_LOWER)
    lower = _invoke(calls, "fixed_y_lower_y44", builder, *lower_args, **{key: value for key, value in lower_kwargs.items() if not key.startswith("_")})
    _require(lower_kwargs["_mu_state"]["calls"] == 1, "mu gate was not called exactly once")
    _require(lower_kwargs["_q_state"]["calls"] == 1, "Q denominator gate was not called exactly once")
    _assert_exact_rows(_rows(_field(lower, "amplitude"), "amplitude", 4), expected_lower["amplitude"], "amplitude")
    _assert_exact_rows(_rows(_field(lower, "phi"), "phi", 4), expected_lower["phi"], "phi")
    _assert_exact_rows(_rows(_field(lower, "scaled_m"), "scaled_m", 4), expected_lower["scaled_m"], "scaled_m")
    _assert_exact_rows(_rows(_field(lower, "scaled_Q"), "scaled_Q", 3), expected_lower["scaled_Q"], "scaled_Q")
    _assert_exact_rows(tuple(_field(lower, "scaled_Qy")), expected_lower["scaled_Qy"], "scaled_Qy")
    _assert_exact_rows(_rows(_field(lower, "reduced_history"), "reduced_history", 4), expected_lower["reduced_history"], "reduced_history")
    _assert_exact_rows(_rows(_field(lower, "source_row_target"), "source_row_target", 4), expected_lower["source_row_target"], "source_row_target")
    _assert_exact_rows(_rows(_field(lower, "incoming_product"), "incoming_product", 4), expected_lower["incoming_product"], "incoming_product")
    _assert_exact_rows(_rows(_field(lower, "selected_product"), "selected_product", 4), expected_lower["selected_product"], "selected_product")
    _require(_field(lower, "fixed_y") == Y_LOWER and _field(lower, "coordinate") == "fixed_y", "lower coordinate changed")
    _require(_field(lower, "log_scale") == expected_lower["log_scale"], "lower log scale changed")
    _require(_field(lower, "log_RQ_scale") == expected_lower["log_RQ_scale"], "lower RQ log scale changed")
    _require(_field(lower, "scaled_RQ_y") == expected_lower["scaled_RQ_y"], "lower RQ mantissa changed")
    _require(_field(lower, "scaled_RQ_y_Z") == expected_lower["scaled_RQ_y_Z"], "lower RQ Z1 changed")
    _assert_uncertified(lower, "lower result")
    _require(_field(lower, "incoming_history_retained") is True, "incoming m1 was dropped")
    _require(_field(lower, "incoming_nonzero_retained") is True, "nonzero incoming m1 was not retained")
    _require(_field(lower, "pstar_applied_once") is True, "Pstar once flag changed")
    _require(_field(lower, "scale_materialized") is False, "common exponential was materialized")

    upper_args, upper_kwargs = _builder_kwargs(Y_UPPER)
    upper = _invoke(calls, "fixed_y_upper_y48", builder, *upper_args, **{key: value for key, value in upper_kwargs.items() if not key.startswith("_")})
    _require(upper_kwargs["_mu_state"]["calls"] == 1 and upper_kwargs["_q_state"]["calls"] == 1, "upper gates were not called exactly once")
    _assert_exact_rows(_rows(_field(upper, "scaled_m"), "upper scaled_m", 4), expected_upper["scaled_m"], "upper scaled_m")
    _assert_exact_rows(_rows(_field(upper, "scaled_Q"), "upper scaled_Q", 3), expected_upper["scaled_Q"], "upper scaled_Q")
    _require(_field(upper, "fixed_y") == Y_UPPER, "upper coordinate changed")
    _require(_field(lower, "scaled_Q") .rows == _field(upper, "scaled_Q").rows, "same-source Q rows differ between endpoints")
    _require(_field(lower, "log_RQ_scale") == _field(upper, "log_RQ_scale"), "RQ/Rp endpoint log scales differ")
    _require(_field(lower, "scaled_RQ_y") == _field(upper, "scaled_RQ_y"), "RQ/Rp endpoint mantissas differ")
    _require(_field(lower, "scaled_RQ_y_Z") == _field(upper, "scaled_RQ_y_Z"), "RQ/Rp endpoint Z1 rows differ")
    _assert_uncertified(upper, "upper result")

    named_mu_gate, named_mu_state = _gate(True)
    named_q_gate, named_q_state = _gate(True)
    named = _invoke(
        calls,
        "named_coordinate_chain_y44",
        builder,
        PSTAR,
        U,
        MU,
        INCOMING_FACTOR1,
        INCOMING_M1,
        AP,
        SCALED_MAIN_ROWS1,
        LOG_END_SCALE,
        mu_nonzero_gate=named_mu_gate,
        q_denominator_nonzero_gate=named_q_gate,
        Z=Z_VALUE,
        delta0=DELTA[0],
        delta1=DELTA[1],
        delta2=DELTA[2],
        fixed_y=Y_LOWER,
    )
    _require(named_mu_state["calls"] == 1 and named_q_state["calls"] == 1, "named coordinate gates changed")
    _require(_rows(_field(named, "scaled_Q"), "named scaled_Q", 3) == expected_lower["scaled_Q"], "named coordinate chain changed Q")
    _require(_field(named, "fixed_y") == Y_LOWER, "named fixed-y chain changed")

    changed_m1 = tuple(value * 2 for value in INCOMING_M1)
    changed_args, changed_kwargs = _builder_kwargs(Y_LOWER, incoming_m1=changed_m1)
    changed = _invoke(calls, "incoming_m1_doubled_changes_result", builder, *changed_args, **{key: value for key, value in changed_kwargs.items() if not key.startswith("_")})
    changed_expected = _expected_core(Y_LOWER, incoming_m1=changed_m1)
    _assert_exact_rows(_rows(_field(changed, "scaled_Q"), "changed-m1 scaled_Q", 3), changed_expected["scaled_Q"], "changed-m1 scaled_Q")
    _require(_field(changed, "incoming_product") .rows != _field(lower, "incoming_product").rows, "changing incoming m1 did not change incoming product")
    _require(_field(changed, "scaled_Q0") != _field(lower, "scaled_Q0"), "changing incoming m1 did not change Q")

    doubled_pstar = tuple(value * 2 for value in PSTAR)
    pstar_args, pstar_kwargs = _builder_kwargs(Y_LOWER, pstar=doubled_pstar)
    pstar_result = _invoke(calls, "Pstar_doubled_once_check", builder, *pstar_args, **{key: value for key, value in pstar_kwargs.items() if not key.startswith("_")})
    pstar_expected = _expected_core(Y_LOWER, pstar=doubled_pstar)
    _assert_exact_rows(_rows(_field(pstar_result, "amplitude"), "doubled-Pstar amplitude", 4), pstar_expected["amplitude"], "doubled-Pstar amplitude")
    _assert_exact_rows(_rows(_field(pstar_result, "scaled_m"), "doubled-Pstar scaled_m", 4), pstar_expected["scaled_m"], "doubled-Pstar scaled_m")
    _assert_exact_rows(_rows(_field(pstar_result, "scaled_Q"), "doubled-Pstar scaled_Q", 3), pstar_expected["scaled_Q"], "doubled-Pstar scaled_Q")
    _require(
        all(left == 2 * right for left, right in zip(_rows(_field(pstar_result, "scaled_Q"), "Pstar Q", 3), expected_lower["scaled_Q"])),
        "Pstar was not applied exactly once",
    )

    zero_mu = (Fraction(0), Fraction(0), Fraction(0), Fraction(0))
    zero_mu_gate, zero_mu_state = _gate(False)
    zero_mu_q_gate, zero_mu_q_state = _gate(True)
    _expected_failure(
        calls,
        "mu_zero_gate_rejected_before_reciprocal",
        builder,
        (ValueError,),
        *(_builder_kwargs(Y_LOWER, mu=zero_mu)[0]),
        mu_nonzero_gate=zero_mu_gate,
        q_denominator_nonzero_gate=zero_mu_q_gate,
    )
    _require(zero_mu_state["calls"] == 1 and zero_mu_q_state["calls"] == 0, "mu=0 gate order changed")

    q_zero_gate, q_zero_state = _gate(False)
    q_zero_mu_gate, q_zero_mu_state = _gate(True)
    q_zero_args = list(lower_args)
    q_zero_args[9] = Fraction(4)
    _expected_failure(
        calls,
        "Q_denominator_zero_gate_rejected",
        builder,
        (ValueError,),
        *tuple(q_zero_args),
        mu_nonzero_gate=q_zero_mu_gate,
        q_denominator_nonzero_gate=q_zero_gate,
    )
    _require(q_zero_mu_state["calls"] == 1 and q_zero_state["calls"] == 1, "Q denominator gate order changed")

    float_args = list(lower_args)
    float_pstar = list(PSTAR)
    float_pstar[0] = 1.0
    float_args[0] = tuple(float_pstar)
    float_mu_gate, float_mu_state = _gate(True)
    float_q_gate, float_q_state = _gate(True)
    _expected_failure(
        calls,
        "builtin_float_rejected",
        builder,
        (TypeError,),
        *tuple(float_args),
        mu_nonzero_gate=float_mu_gate,
        q_denominator_nonzero_gate=float_q_gate,
    )
    _require(float_mu_state["calls"] == 0 and float_q_state["calls"] == 0, "float was not rejected before gates")

    symbolic_t = SymbolicScalar.symbol("g_gap_one_minus_exp_minus_4_symbolic")
    q_rows = expected_lower["scaled_Q"]
    k_rows = expected_lower["K"]
    functional_expected = _expected_functional(q_rows, k_rows, symbolic_t)
    functional_q_gate, functional_q_state = _gate(True)
    functional = _invoke(
        calls,
        "conditional_gap_functional_symbolic_decay_Z1",
        signed,
        q_rows,
        Z_VALUE,
        DELTA[0],
        DELTA[1],
        k_rows,
        Fraction(0),
        Fraction(0),
        Y_LOWER,
        Y_UPPER,
        Fraction(0),
        Fraction(0),
        symbolic_t,
        q_denominator_nonzero_gate=functional_q_gate,
    )
    _require(functional_q_state["calls"] == 1, "functional Q gate was not called exactly once")
    for name in (
        "linear_log_scale",
        "linear_log_scale_Z",
        "linear_mantissa",
        "linear_physical_Z_mantissa",
        "quadratic_log_scale",
        "quadratic_log_scale_Z",
        "quadratic_mantissa",
        "quadratic_physical_Z_mantissa",
        "c",
        "c_Z",
        "gap_decay_Z",
    ):
        _require(_field(functional, name) == functional_expected[name], "functional hand check changed: " + name)
    _require(_field(functional, "gap_decay") == symbolic_t, "gap decay was materialized")
    _assert_uncertified(functional, "signed functional")

    checks = {
        "nondegenerate_gap_xi_11_to_12": XI_LOWER < XI_UPPER and (MU[0] * Y_LOWER, MU[0] * Y_UPPER) == (XI_LOWER, XI_UPPER),
        "fixed_y_endpoint_pair_y44_y48": _field(lower, "fixed_y") == Y_LOWER and _field(upper, "fixed_y") == Y_UPPER,
        "ordinary_Z0_to_Z3_exact": len(_rows(_field(lower, "scaled_m"), "final scaled_m", 4)) == 4,
        "same_source_RQy_Rp_endpoint": _field(lower, "scaled_RQ_y") == _field(upper, "scaled_RQ_y"),
        "same_source_RQy_Rp_Z1": _field(lower, "scaled_RQ_y_Z") == _field(upper, "scaled_RQ_y_Z"),
        "incoming_m1_nonzero_retained": all(value != 0 for value in INCOMING_M1),
        "incoming_m1_change_changes_result": _field(changed, "scaled_Q0") != _field(lower, "scaled_Q0"),
        "Pstar_applied_exactly_once": _field(lower, "pstar_applied_once") is True and all(left == 2 * right for left, right in zip(_rows(_field(pstar_result, "scaled_Q"), "Pstar Q final", 3), expected_lower["scaled_Q"])),
        "mu_zero_gate": calls[5]["passed"],
        "Q_denominator_zero_gate": calls[6]["passed"],
        "builtin_float_rejected": calls[7]["passed"],
        "coordinate_chain_named_equals_positional": _rows(_field(named, "scaled_Q"), "named final Q", 3) == expected_lower["scaled_Q"],
        "symbolic_gap_decay_not_materialized": isinstance(_field(functional, "gap_decay"), SymbolicScalar),
        "ordinary_Z1_functional_hand_check": _field(functional, "linear_physical_Z_mantissa") == functional_expected["linear_physical_Z_mantissa"] and _field(functional, "quadratic_physical_Z_mantissa") == functional_expected["quadratic_physical_Z_mantissa"],
        "source_flags_not_certified": all(_field(lower, key) is False for key in ("source_identity_proved", "physical_scale_proved", "certified", "error_certified")),
    }
    _require(all(checks.values()), "reduced-history fixture contains a failed check")
    return {
        "schema": COMPONENT + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_checks_passed",
        "passed": True,
        "checks_passed": True,
        "synthetic_only": True,
        "call_count": len(calls),
        "calls": calls,
        "ordinary_Z_orders": [0, 1, 2, 3],
        "gap": {
            "mu": MU[0],
            "xi_lower": XI_LOWER,
            "xi_upper": XI_UPPER,
            "y_lower": Y_LOWER,
            "y_upper": Y_UPPER,
            "nondegenerate": True,
            "endpoints_are_synthetic_coordinates": True,
            "physical_source_value_claim": False,
        },
        "hand_calculation": {
            "amplitude_rows_y44": expected_lower["amplitude"],
            "scaled_Q_rows_y44": expected_lower["scaled_Q"],
            "scaled_Qy_rows_y44": expected_lower["scaled_Qy"],
            "functional_c": functional_expected["c"],
            "functional_c_Z": functional_expected["c_Z"],
            "functional_linear_log_scale": functional_expected["linear_log_scale"],
            "functional_linear_mantissa": functional_expected["linear_mantissa"],
            "functional_linear_physical_Z_mantissa": functional_expected["linear_physical_Z_mantissa"],
            "functional_quadratic_log_scale": functional_expected["quadratic_log_scale"],
            "functional_quadratic_mantissa": functional_expected["quadratic_mantissa"],
            "functional_quadratic_physical_Z_mantissa": functional_expected["quadratic_physical_Z_mantissa"],
            "gap_decay_symbolic_only": True,
        },
        "expected_lower": expected_lower,
        "expected_upper": expected_upper,
        "expected_integral": functional_expected,
        "integral_decay": {
            "symbol": symbolic_t,
            "meaning": "caller-supplied g=1-exp(-4), symbolic only; not a materialized value",
            "materialized": False,
        },
        "same_source_endpoint_premises": {
            "same_caller_rows": True,
            "same_mu_delta_Pstar": True,
            "Rp_Z_zero_premise": True,
            "Z1_comparison_is_conditional_algebra": True,
            "source_identity_proved": False,
        },
        "checks": checks,
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "physical_scale_proved": False,
        "physical_scale_certified": False,
        "certification": False,
        "certified": False,
        "error_certified": False,
        "error_terms_deletable": False,
        "uncertified": True,
    }


def _expected_functional(q: Jet3, k: Jet4, t: SymbolicScalar) -> dict[str, Any]:
    z = Z_VALUE
    d0, d1 = DELTA[0], DELTA[1]
    L = Fraction(1) - d0 * z * z
    L_Z = -(d1 * z * z + d0 * 2 * z)
    numerator = (Fraction(1) - d0) * z
    numerator_Z = -d1 * z + (Fraction(1) - d0)
    c = numerator / (Fraction(4) * L)
    c_Z = (numerator_Z * L - numerator * L_Z) / (Fraction(4) * L * L)
    q0, q1, q2 = q
    linear_mantissa = c * q1 * t
    linear_physical_Z_mantissa = c_Z * q1 * t + c * (q2 + q1 * Fraction(0)) * t
    quadratic_mantissa = -(q0 * q0)  # delta-y = 48-44 = 4; /4 cancels it.
    quadratic_physical_Z_mantissa = -(2 * q0 * q1)
    return {
        "linear_log_scale": k[0] + Y_UPPER,
        "linear_log_scale_Z": k[1],
        "linear_mantissa": linear_mantissa,
        "linear_physical_Z_mantissa": linear_physical_Z_mantissa,
        "quadratic_log_scale": 2 * k[0],
        "quadratic_log_scale_Z": 2 * k[1],
        "quadratic_mantissa": quadratic_mantissa,
        "quadratic_physical_Z_mantissa": quadratic_physical_Z_mantissa,
        "c": c,
        "c_Z": c_Z,
        "gap_decay_Z": SymbolicScalar(0),
    }


def _walk_forbidden_true(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                found.append(path + "." + key_text)
            found.extend(_walk_forbidden_true(item, path + "." + key_text))
    elif isinstance(value, (list, tuple, set, frozenset)):
        for index, item in enumerate(value):
            found.extend(_walk_forbidden_true(item, path + "[" + str(index) + "]"))
    return found


def _contains_builtin_float(value: Any) -> bool:
    if type(value) is float:
        return True
    if isinstance(value, Mapping):
        return any(_contains_builtin_float(item) for item in value.values())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_builtin_float(item) for item in value)
    if isinstance(value, SymbolicScalar):
        return any(_contains_builtin_float(item) for _, item in value._terms)
    if isinstance(value, Fraction):
        return False
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(_contains_builtin_float(getattr(value, field.name)) for field in dataclasses.fields(value))
    fields = getattr(value, "__dict__", None)
    return isinstance(fields, dict) and _contains_builtin_float(fields)


def _contains_midpoint(value: Any) -> bool:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if "midpoint" in str(key).lower():
                return True
            if _contains_midpoint(item):
                return True
        return False
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_midpoint(item) for item in value)
    if isinstance(value, SymbolicScalar) or isinstance(value, Fraction):
        return False
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return any(_contains_midpoint(getattr(value, field.name)) for field in dataclasses.fields(value))
    fields = getattr(value, "__dict__", None)
    return isinstance(fields, dict) and _contains_midpoint(fields)


def _validate_fixture(value: Mapping[str, Any]) -> None:
    _require(value.get("passed") is True and value.get("checks_passed") is True, "fixture did not pass")
    _require(value.get("synthetic_only") is True, "fixture is not synthetic-only")
    calls = value.get("calls")
    _require(isinstance(calls, list) and len(calls) == value.get("call_count"), "fixture call log is incomplete")
    _require(all(isinstance(item, Mapping) and item.get("invoked") is True and item.get("passed") is True for item in calls), "fixture call log contains a failed call")
    _require(not _walk_forbidden_true(value), "fixture made a forbidden source or physical claim")
    _require(not _contains_builtin_float(value), "fixture exposed a built-in float")
    _require(not _contains_midpoint(value), "fixture exposed an interval midpoint")
    checks = value.get("checks")
    _require(isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()), "fixture contains a failed check")


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
    preflight: Mapping[str, Any],
    leaf_imported: bool,
    error: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    flags = {
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "source_values_computed": False,
        "physical_scale_proved": False,
        "physical_scale_certified": False,
        "physical_source_owner_proved": False,
        "source_Rp_owner_authenticated": False,
        "source_parameter_independence_proved": False,
        "source_mu_independence_proved": False,
        "source_delta_independence_proved": False,
        "source_Rp_independence_proved": False,
        "error_certified": False,
        "error_terms_deletable": False,
        "numeric_use_admissible": False,
        "certification": False,
        "certified": False,
        "old_math_called": False,
        "source_graph_traversed": False,
        "ancestor_replayed": False,
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
        "preflight": dict(preflight),
        "raw_sha256_authenticated_before_entry": True,
        "ast_parse_before_import": True,
        "leaf_imported": leaf_imported,
        "fixture_invocations": 1 if leaf_imported else 0,
        "leaf_api_invocations": len(calls),
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
            "same_source_endpoint_comparison_conditional": True,
            "source_owned_F_bound": False,
            "physical_conclusions": False,
            "uncertified": True,
        },
        "error": _json_safe(error),
    }


def main() -> int:
    source, hashes, source_sha, brief_sha, preflight = _preflight()
    ast_info = _static_parse_leaf()
    entry_sha = _write_exclusive(
        ENTERED,
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
            "leaf_api_invocations": 0,
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
            preflight=preflight,
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
            preflight=preflight,
            leaf_imported=leaf_imported,
            error={
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        )

    result_sha = _write_atomic(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "entry_sha256": entry_sha,
                "fixture_invocations": output["fixture_invocations"],
                "leaf_api_invocations": output["leaf_api_invocations"],
                "source_identity_proved": False,
                "physical_scale_proved": False,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
