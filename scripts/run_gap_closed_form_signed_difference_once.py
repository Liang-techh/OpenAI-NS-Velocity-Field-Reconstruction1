"""Run the closed-form signed-gap difference leaf exactly once.

This runner authenticates only the new leaf and this runner by raw SHA-256.
After an exclusive ``entered.json`` receipt is written, it imports and calls
only the new leaf.  The fixture is synthetic exact arithmetic: it checks the
four independently scaled rows, the complete physical-Z1 ports, the
correlated gap-decay identities, and fail-closed input guards.  It does not
read a source archive, call an ancestor, traverse a graph, or claim source,
physical-scale, interval, quadrature, or error certification.
"""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "gap_closed_form_signed_difference_v1"
PHASE = COMPONENT + "_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/gap_closed_form_signed_difference.py"
RUNNER_REL = "scripts/run_gap_closed_form_signed_difference_once.py"
EXPECTED_FILES = frozenset((LEAF_REL, RUNNER_REL))
LEAF = ROOT / LEAF_REL
RUNNER = ROOT / RUNNER_REL

F = Fraction


class RunnerInputError(RuntimeError):
    """Raised when the pending one-entry contract is not frozen."""


class MidpointFixture:
    """An explicit midpoint marker used only by the negative control."""

    is_midpoint = True

    def __repr__(self) -> str:
        return "MidpointFixture()"


@dataclasses.dataclass(frozen=True, slots=True)
class DirectedIntervalFixture:
    """Exact directed bounds used by the single interval consumer case."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction) or not isinstance(self.upper, Fraction):
            raise TypeError("interval bounds must be Fractions")
        if self.lower > self.upper:
            raise ValueError("directed interval bounds are reversed")

    def __mul__(self, other: Any) -> Any:
        if isinstance(other, Fraction):
            if other >= 0:
                return DirectedIntervalFixture(self.lower * other, self.upper * other)
            return DirectedIntervalFixture(self.upper * other, self.lower * other)
        if isinstance(other, DirectedIntervalFixture):
            products = (
                self.lower * other.lower,
                self.lower * other.upper,
                self.upper * other.lower,
                self.upper * other.upper,
            )
            return DirectedIntervalFixture(min(products), max(products))
        return NotImplemented

    def __rmul__(self, other: Any) -> Any:
        return self * other

    def __repr__(self) -> str:
        return f"DirectedIntervalFixture([{self.lower},{self.upper}])"


SIGNED_ROWS: dict[str, dict[str, Fraction]] = {
    "linear": {
        "LF": F(5, 2),
        "aF": F(3),
        "delta_a": F(2),
        "aF_physZ": F(7),
        "delta_physZ": F(5),
        "t": F(1, 2),
    },
    "quadratic": {
        "LF": F(-7, 3),
        "aF": F(-4),
        "delta_a": F(3),
        "aF_physZ": F(2),
        "delta_physZ": F(-5),
        "t": F(-1, 3),
    },
    "endpoint_upper": {
        "LF": F(11, 4),
        "aF": F(6),
        "delta_a": F(-2),
        "aF_physZ": F(-3),
        "delta_physZ": F(4),
        "t": F(3, 4),
    },
    "endpoint_lower": {
        "LF": F(-13, 5),
        "aF": F(-5),
        "delta_a": F(4),
        "aF_physZ": F(8),
        "delta_physZ": F(-6),
        "t": F(2, 5),
    },
}

SIGNED_EXPECTED: dict[str, dict[str, Any]] = {
    "linear": {
        "label": "linear+",
        "sign": "+",
        "LF": F(5, 2),
        "value": F(9, 2),
        "physicalZ": F(11),
        "contribution": F(9, 2),
        "contribution_physicalZ": F(11),
    },
    "quadratic": {
        "label": "quadratic+",
        "sign": "+",
        "LF": F(-7, 3),
        "value": F(10, 3),
        "physicalZ": F(-4),
        "contribution": F(10, 3),
        "contribution_physicalZ": F(-4),
    },
    "endpoint_upper": {
        "label": "upper-",
        "sign": "-",
        "LF": F(11, 4),
        "value": F(1),
        "physicalZ": F(19, 4),
        "contribution": F(-1),
        "contribution_physicalZ": F(-19, 4),
    },
    "endpoint_lower": {
        "label": "lower+",
        "sign": "+",
        "LF": F(-13, 5),
        "value": F(18, 5),
        "physicalZ": F(-26, 5),
        "contribution": F(18, 5),
        "contribution_physicalZ": F(-26, 5),
    },
}


GAP_ROWS: dict[str, dict[str, Fraction]] = {
    "linear": {"hF": F(3, 5), "d": F(2, 7), "etaF": F(1, 3), "deta": F(5, 11)},
    "quadratic": {"hF": F(4, 7), "d": F(-1, 5), "etaF": F(2, 9), "deta": F(3, 10)},
    "endpoint_upper": {"hF": F(5, 6), "d": F(3, 8), "etaF": F(-1, 4), "deta": F(7, 13)},
    "endpoint_lower": {"hF": F(7, 8), "d": F(-2, 9), "etaF": F(1, 5), "deta": F(-4, 15)},
}

GAP_EXPECTED: dict[str, dict[str, Any]] = {
    "linear": {
        "label": "linear+",
        "sign": "+",
        "delta_g": F(6, 35),
        "delta_g_Z": F(53, 385),
        "contribution": F(6, 35),
        "contribution_Z": F(53, 385),
    },
    "quadratic": {
        "label": "quadratic+",
        "sign": "+",
        "delta_g": F(-4, 35),
        "delta_g_Z": F(364, 1575),
        "contribution": F(-4, 35),
        "contribution_Z": F(364, 1575),
    },
    "endpoint_upper": {
        "label": "upper-",
        "sign": "-",
        "delta_g": F(5, 16),
        "delta_g_Z": F(895, 2496),
        "contribution": F(-5, 16),
        "contribution_Z": F(-895, 2496),
    },
    "endpoint_lower": {
        "label": "lower+",
        "sign": "+",
        "delta_g": F(-7, 36),
        "delta_g_Z": F(-133, 540),
        "contribution": F(-7, 36),
        "contribution_Z": F(-133, 540),
    },
}


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
        return str(value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _json_safe(getattr(value, field.name))
            for field in dataclasses.fields(value)
        }
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=repr)
    if isinstance(value, Path):
        return value.as_posix()
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


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_bytes())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _repo_path(raw_path: Any, label: str) -> tuple[str, Path]:
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
    files = source.get("files")
    _require(isinstance(files, list) and len(files) == 2, "source manifest must contain exactly two files")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_path(item.get("path"), f"source.files[{index}]")
        _require(key in EXPECTED_FILES, "source manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "source manifest contains duplicate file: " + key)
        expected = item.get("sha256")
        _require(isinstance(expected, str) and len(expected) == 64, "invalid SHA-256 for " + key)
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        hashes[key] = actual
    _require(frozenset(hashes) == EXPECTED_FILES, "manifest must contain exactly leaf and runner")
    return hashes


def _brief_scope(raw: Mapping[str, Any]) -> dict[str, Any]:
    if raw.get("component") == COMPONENT:
        return dict(raw)
    for key in ("pending_components", "components", "pending"):
        table = raw.get(key)
        if isinstance(table, Mapping) and isinstance(table.get(COMPONENT), Mapping):
            return dict(table[COMPONENT])
    raise RunnerInputError("real pending brief has no component section: " + COMPONENT)


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(DIRECTORY.is_dir(), "component directory is missing")
    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(BRIEF.is_file(), "real pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")

    brief = _brief_scope(_load_json_object(BRIEF, "real pending brief"))
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal guards changed",
    )
    source_path = SOURCE.resolve()
    checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    checkpoint_path = Path(str(checkpoint.get("path", "")).replace("\\", "/"))
    if not checkpoint_path.is_absolute():
        checkpoint_path = ROOT / checkpoint_path
    source_sha = _sha256_file(source_path)
    _require(checkpoint_path.resolve() == source_path, "pending brief checkpoint path mismatch")
    _require(
        isinstance(checkpoint.get("sha256"), str)
        and checkpoint["sha256"].lower() == source_sha,
        "pending brief checkpoint SHA mismatch",
    )

    source = _load_json_object(SOURCE, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("schema") == COMPONENT + "_checkpoint_v1", "source manifest schema mismatch")
    _require(source.get("one_entry_only") is True, "one-entry guard missing")
    _require(source.get("synthetic_only") is True, "synthetic-only guard missing")
    for key in ("source_identity_proved", "physical_scale_proved", "error_certified"):
        _require(source.get(key) is False, key + " must remain false")
    if "certified" in source:
        _require(source["certified"] is False, "certified must remain false")
    if "no_ancestor_replay" in source:
        _require(source["no_ancestor_replay"] is True, "ancestor replay guard changed")
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha256_file(BRIEF)


def _import_leaf() -> Any:
    _require(LEAF.is_file(), "new signed-gap leaf is missing")
    spec = importlib.util.spec_from_file_location("gap_closed_form_signed_difference_v1_new_leaf", LEAF)
    _require(spec is not None and spec.loader is not None, "new signed-gap leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "result is missing field: " + name)
        return value[name]
    _require(hasattr(value, name), "result is missing field: " + name)
    return getattr(value, name)


def _interval_bounds(value: Any, label: str) -> tuple[Fraction, Fraction]:
    _require(isinstance(value, DirectedIntervalFixture), label + " is not a directed interval")
    return value.lower, value.upper


def _call(
    calls: list[dict[str, Any]],
    *,
    name: str,
    function_name: str,
    function: Callable[..., Any],
    expected_exception: type[BaseException] | None = None,
    **kwargs: Any,
) -> Any:
    record: dict[str, Any] = {
        "name": name,
        "function": function_name,
        "kwargs": _json_safe(kwargs),
        "invoked": False,
        "expected_exception": expected_exception.__name__ if expected_exception else None,
    }
    calls.append(record)
    try:
        record["invoked"] = True
        output = function(**kwargs)
        record["returned"] = True
        record["output"] = _json_safe(output)
        if expected_exception is not None:
            record["passed"] = False
            raise RunnerInputError(name + " unexpectedly returned")
        record["passed"] = True
        return output
    except BaseException as exc:
        record["returned"] = False
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["traceback"] = traceback.format_exc()
        if expected_exception is None:
            record["passed"] = False
            raise
        record["passed"] = isinstance(exc, expected_exception)
        if not record["passed"]:
            raise
        return None


def _run_signed_fixture(leaf: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    builder = getattr(leaf, "build_gap_closed_form_signed_difference", None)
    _require(callable(builder), "signed closed-form builder is missing")
    error_type = getattr(leaf, "GapClosedFormSignedDifferenceInputError", None)
    _require(isinstance(error_type, type), "signed closed-form input error is missing")

    result = _call(
        calls,
        name="four_independent_signed_rows",
        function_name="build_gap_closed_form_signed_difference",
        function=builder,
        rows=SIGNED_ROWS,
    )
    terms = _field(result, "terms")
    _require(isinstance(terms, tuple) and len(terms) == 4, "four signed rows were not retained")
    checks: dict[str, bool] = {}
    for row_name, expected in SIGNED_EXPECTED.items():
        component = _field(result, row_name)
        checks[row_name + "_label"] = _field(component, "name") == expected["label"]
        checks[row_name + "_sign"] = _field(component, "sign") == expected["sign"]
        checks[row_name + "_scale_retained"] = _field(component, "LF") == expected["LF"]
        checks[row_name + "_value"] = _field(component, "signed_value_mantissa") == expected["value"]
        checks[row_name + "_physZ"] = _field(component, "signed_physicalZ_mantissa") == expected["physicalZ"]
        checks[row_name + "_contribution"] = _field(component, "contribution_mantissa") == expected["contribution"]
        checks[row_name + "_contribution_physZ"] = _field(component, "contribution_physicalZ_mantissa") == expected["contribution_physicalZ"]
        checks[row_name + "_direct_delta_nonzero"] = SIGNED_ROWS[row_name]["delta_a"] != 0
        checks[row_name + "_t_nonzero"] = SIGNED_ROWS[row_name]["t"] != 0
        checks[row_name + "_phys_ports_distinct"] = (
            SIGNED_ROWS[row_name]["aF_physZ"] != 0
            and SIGNED_ROWS[row_name]["delta_physZ"] != 0
        )
    checks.update(
        {
            "ordinary_Z1": _field(result, "ordinary_Z1") is True,
            "scales_not_materialized": _field(result, "scales_materialized") is False,
            "endpoints_retained": _field(result, "endpoint_terms_retained") is True,
            "terms_not_combined": _field(result, "terms_not_combined") is True,
            "frame_logscale_not_materialized": _field(result, "frame_logscale_materialized") is False,
            "t_expm1_caller_premise": _field(result, "t_expm1_caller_premise") is True,
            "exponential_not_evaluated": _field(result, "exponential_evaluated") is False,
            "direct_difference_caller_premise": _field(result, "direct_difference_identity_caller_premise") is True,
            "LF_Z_not_added": _field(result, "LF_Z_added") is False,
            "t_Z_not_added": _field(result, "t_Z_added") is False,
            "P1_is_not_m15": _field(result, "P1_is_m15") is False,
            "P1_not_m15": _field(result, "P1_not_m15") is True,
        }
    )
    for name in (
        "source_ownership",
        "F_ownership",
        "source_identity_proved",
        "same_source_proved",
        "physical_scale_proved",
        "interval_certified",
        "certified",
        "quadrature_error_certified",
        "error_certified",
        "direct_difference_identity_proved",
    ):
        checks[name + "_false"] = _field(result, name) is False
    _require(all(checks.values()), "signed closed-form fixture check failed")

    _call(
        calls,
        name="signed_missing_row_rejected",
        function_name="build_gap_closed_form_signed_difference",
        function=builder,
        expected_exception=error_type,
        rows={name: value for name, value in SIGNED_ROWS.items() if name != "endpoint_lower"},
    )
    midpoint_row = dict(SIGNED_ROWS["linear"])
    midpoint_row["LF"] = MidpointFixture()
    _call(
        calls,
        name="signed_midpoint_rejected",
        function_name="build_gap_closed_form_signed_difference",
        function=builder,
        expected_exception=error_type,
        rows={**SIGNED_ROWS, "linear": midpoint_row},
    )
    float_row = dict(SIGNED_ROWS["linear"])
    float_row["LF"] = 1.0
    _call(
        calls,
        name="signed_float_rejected",
        function_name="build_gap_closed_form_signed_difference",
        function=builder,
        expected_exception=error_type,
        rows={**SIGNED_ROWS, "linear": float_row},
    )
    return {"checks": checks, "result": _json_safe(result)}


def _run_gap_fixture(leaf: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    builder = getattr(leaf, "build_correlated_gap_decay_difference", None)
    _require(callable(builder), "correlated gap-decay builder is missing")
    error_type = getattr(leaf, "GapClosedFormSignedDifferenceInputError", None)
    _require(isinstance(error_type, type), "gap-decay input error is missing")

    result = _call(
        calls,
        name="four_correlated_gap_rows",
        function_name="build_correlated_gap_decay_difference",
        function=builder,
        rows=GAP_ROWS,
    )
    terms = _field(result, "terms")
    _require(isinstance(terms, tuple) and len(terms) == 4, "four gap rows were not retained")
    checks: dict[str, bool] = {}
    for row_name, expected in GAP_EXPECTED.items():
        component = _field(result, row_name)
        checks[row_name + "_label"] = _field(component, "name") == expected["label"]
        checks[row_name + "_sign"] = _field(component, "sign") == expected["sign"]
        checks[row_name + "_delta_g"] = _field(component, "signed_gap_mantissa") == expected["delta_g"]
        checks[row_name + "_delta_g_Z"] = _field(component, "signed_gap_physZ_mantissa") == expected["delta_g_Z"]
        checks[row_name + "_contribution"] = _field(component, "contribution_mantissa") == expected["contribution"]
        checks[row_name + "_contribution_Z"] = _field(component, "contribution_physZ_mantissa") == expected["contribution_Z"]
        checks[row_name + "_nonzero_delta_g"] = expected["delta_g"] != 0
        checks[row_name + "_nonzero_delta_g_Z"] = expected["delta_g_Z"] != 0
    checks.update(
        {
            "ordinary_Z1": _field(result, "ordinary_Z1") is True,
            "direct_gap_ports_caller_premise": _field(result, "direct_gap_ports_caller_premise") is True,
            "exponential_not_evaluated": _field(result, "exponential_evaluated") is False,
            "endpoints_retained": _field(result, "endpoint_terms_retained") is True,
            "terms_not_combined": _field(result, "terms_not_combined") is True,
        }
    )
    for name in (
        "source_ownership",
        "F_ownership",
        "source_identity_proved",
        "same_source_proved",
        "physical_scale_proved",
        "interval_certified",
        "certified",
        "quadrature_error_certified",
        "error_certified",
        "direct_gap_ports_proved",
    ):
        checks[name + "_false"] = _field(result, name) is False
    _require(all(checks.values()), "correlated gap-decay fixture check failed")

    interval_rows = dict(GAP_ROWS)
    interval_rows["linear"] = {
        "hF": DirectedIntervalFixture(F(0), F(1, 16)),
        "d": F(1, 2),
        "etaF": F(1),
        "deta": F(1),
    }
    interval_result = _call(
        calls,
        name="directed_interval_gap_row",
        function_name="build_correlated_gap_decay_difference",
        function=builder,
        rows=interval_rows,
    )
    interval_linear = _field(interval_result, "linear")
    interval_delta_g = _interval_bounds(
        _field(interval_linear, "signed_gap_mantissa"),
        "interval Delta g",
    )
    interval_delta_g_z = _interval_bounds(
        _field(interval_linear, "signed_gap_physZ_mantissa"),
        "interval Delta g_Z",
    )
    interval_input_values = tuple(interval_rows["linear"].values())
    interval_checks = {
        "hF_bounds": interval_rows["linear"]["hF"] == DirectedIntervalFixture(F(0), F(1, 16)),
        "delta_g_bounds": interval_delta_g == (F(0), F(1, 32)),
        "delta_g_Z_bounds": interval_delta_g_z == (F(0), F(0)),
        "ordinary_Z1": _field(interval_result, "ordinary_Z1") is True,
        "linear_plus": _field(interval_linear, "sign") == "+",
        "interval_not_certified": _field(interval_result, "interval_certified") is False,
        "no_midpoint_or_float": not any(
            isinstance(value, (MidpointFixture, int, float, complex, bool))
            for value in interval_input_values
        ),
    }
    _require(all(interval_checks.values()), "directed interval gap fixture check failed")

    _call(
        calls,
        name="gap_missing_row_rejected",
        function_name="build_correlated_gap_decay_difference",
        function=builder,
        expected_exception=error_type,
        rows={name: value for name, value in GAP_ROWS.items() if name != "endpoint_lower"},
    )
    midpoint_row = dict(GAP_ROWS["linear"])
    midpoint_row["hF"] = MidpointFixture()
    _call(
        calls,
        name="gap_midpoint_rejected",
        function_name="build_correlated_gap_decay_difference",
        function=builder,
        expected_exception=error_type,
        rows={**GAP_ROWS, "linear": midpoint_row},
    )
    float_row = dict(GAP_ROWS["linear"])
    float_row["hF"] = 1.0
    _call(
        calls,
        name="gap_float_rejected",
        function_name="build_correlated_gap_decay_difference",
        function=builder,
        expected_exception=error_type,
        rows={**GAP_ROWS, "linear": float_row},
    )
    return {
        "checks": checks,
        "interval_checks": interval_checks,
        "result": _json_safe(result),
        "interval_result": _json_safe(interval_result),
    }


def main() -> int:
    source: dict[str, Any] | None = None
    hashes: dict[str, str] = {}
    source_sha: str | None = None
    brief_sha: str | None = None
    entry_sha: str | None = None
    calls: list[dict[str, Any]] = []
    leaf_imported = False
    invocation_count = 0

    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": dict(hashes),
            "source_manifest": _json_safe(source),
            "leaf_sha256": hashes.get(LEAF_REL),
            "runner_sha256": hashes.get(RUNNER_REL),
            "leaf_imported": False,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_computed": False,
            "runtime_model_metadata": None,
        },
    )

    try:
        leaf = _import_leaf()
        leaf_imported = True
        signed = _run_signed_fixture(leaf, calls)
        gap = _run_gap_fixture(leaf, calls)
        invocation_count = len(calls)
        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": dict(hashes),
            "source_manifest": _json_safe(source),
            "leaf_imported": leaf_imported,
            "consumer_calls": invocation_count,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "interval_certified": False,
            "certified": False,
            "error_certified": False,
            "synthetic_only": True,
            "runtime_model_metadata": None,
            "result": {"signed": signed, "gap_decay": gap},
            "consumer_call_log": calls,
            "error": None,
        }
    except BaseException as exc:
        invocation_count = len(calls)
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": dict(hashes),
            "source_manifest": _json_safe(source) if source is not None else None,
            "leaf_sha256": hashes.get(LEAF_REL),
            "runner_sha256": hashes.get(RUNNER_REL),
            "leaf_imported": leaf_imported,
            "consumer_calls": invocation_count,
            "ancestor_math_calls": 0,
            "source_graph_traversed": False,
            "source_values_selected": False,
            "source_values_computed": False,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "interval_certified": False,
            "certified": False,
            "error_certified": False,
            "synthetic_only": True,
            "runtime_model_metadata": None,
            "result": None,
            "consumer_call_log": calls,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        }

    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "source_sha256": source_sha,
                "entry_sha256": entry_sha,
                "consumer_calls": invocation_count,
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
