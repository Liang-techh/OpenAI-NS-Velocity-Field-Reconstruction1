"""Run the synthetic F-gap first-row reduction fixture exactly once.

The manifest authenticates only this runner and the new reduction leaf.  The
fixture uses exact rational ordinary-Z rows (orders 0 through 3), performs an
independent product calculation, and keeps source identity and numerical
certification explicitly false.  No producer, graph, ancestor, or old helper
is imported.
"""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.util
import json
import os
import re
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_first_row_reduction_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/f_gap_selected_first_row_reduction.py"
RUNNER_REL = "scripts/run_f_gap_first_row_reduction_once.py"
LEAF = ROOT / LEAF_REL
RUNNER = ROOT / RUNNER_REL
REQUIRED_INPUTS = frozenset({LEAF_REL, RUNNER_REL})
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
ZERO = (Fraction(0), Fraction(0), Fraction(0), Fraction(0))

FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "source_identity_proved",
        "source_values_authenticated",
        "determinant_positivity_proved",
        "error_certified",
        "numeric_certification",
        "numeric_certified",
        "certified",
        "source_graph_traversed",
        "old_source_graph_traversed",
        "ancestor_replayed",
        "old_math_called",
        "physical_conclusion_proved",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


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
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if type(value) is float:
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
    """Write a complete result packet, then publish it with one rename."""

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
    files = source.get("files")
    _require(
        isinstance(files, list) and len(files) == 2,
        "source manifest must contain exactly two file entries",
    )
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_path(item.get("path"), f"source.files[{index}]")
        _require(key in REQUIRED_INPUTS, "manifest contains an out-of-scope file: " + key)
        _require(key not in hashes, "duplicate manifest path: " + key)
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


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path, expected = checkpoint.get("path"), checkpoint.get("sha256")
    _require(
        isinstance(raw_path, str) and isinstance(expected, str),
        "pending brief source checkpoint is incomplete",
    )
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(
        candidate.resolve() == SOURCE.resolve(),
        "pending brief points to a different source checkpoint",
    )
    actual = _sha256_file(SOURCE)
    _require(
        SHA256_RE.fullmatch(expected) is not None and actual.lower() == expected.lower(),
        "pending brief source checkpoint SHA-256 mismatch",
    )
    return actual


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTERED.exists() and not RESULT.exists(), "one-entry component already entered")
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
    _require(
        source.get("schema") == COMPONENT + "_checkpoint_v1",
        "source manifest schema mismatch",
    )
    for key in (
        "one_entry_only",
        "synthetic_only",
        "no_source_graph_traversal",
        "no_ancestor_replay",
    ):
        _require(source.get(key) is True, "source manifest guard missing or changed: " + key)
    for key in FORBIDDEN_TRUE_KEYS:
        if key in source:
            _require(source[key] is False, "source manifest made forbidden claim: " + key)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, hashlib.sha256(brief_raw).hexdigest()


def _import_leaf() -> Any:
    _require(LEAF.is_file(), "new first-row reduction leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_selected_first_row_reduction_new_leaf", LEAF
    )
    _require(spec is not None and spec.loader is not None, "new leaf load failed")
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
    record: dict[str, Any] = {
        "name": name,
        "invoked": True,
        "returned": False,
        "passed": False,
    }
    calls.append(record)
    try:
        value = function(*args, **kwargs)
    except BaseException as exc:
        record["error"] = {"type": type(exc).__name__, "message": str(exc)}
        record["traceback"] = traceback.format_exc()
        raise
    record["returned"] = True
    record["passed"] = True
    record["output"] = _json_safe(value)
    return value


def _rows(value: Any, label: str) -> tuple[Any, ...]:
    raw = value.get("rows") if isinstance(value, Mapping) else getattr(value, "rows", None)
    _require(raw is not None, label + " has no ordinary rows")
    result = tuple(raw)
    _require(len(result) == 4, label + " does not have rows 0 through 3")
    return result


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "leaf result omitted " + name)
        return value[name]
    result = getattr(value, name, None)
    _require(result is not None, "leaf result omitted " + name)
    return result


def _ordinary_product(
    left: tuple[Fraction, Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Independent ordinary derivative product, including explicit Z3."""

    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
        left[3] * right[0]
        + 3 * left[2] * right[1]
        + 3 * left[1] * right[2]
        + left[0] * right[3],
    )


def _ordinary_add(
    left: tuple[Fraction, Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return tuple(left[index] + right[index] for index in range(4))  # type: ignore[return-value]


def _ordinary_neg(
    value: tuple[Fraction, Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return tuple(-item for item in value)  # type: ignore[return-value]


def _rows_are_fractional(rows: tuple[Any, ...], label: str) -> None:
    _require(len(rows) == 4, label + " must have four rows")
    _require(all(type(item) is Fraction for item in rows), label + " is not exact rational data")


def _assert_uncertified(value: Any, label: str) -> None:
    for name in (
        "source_identity_proved",
        "source_values_authenticated",
        "determinant_positivity_proved",
        "error_certified",
    ):
        _require(_field(value, name) is False, label + " made a forbidden claim: " + name)


HAND_INPUTS: dict[str, tuple[Fraction, Fraction, Fraction, Fraction]] = {
    "incoming_factor1": (Fraction(1), Fraction(2), Fraction(3), Fraction(4)),
    "incoming_m1": (Fraction(2), Fraction(3), Fraction(5), Fraction(7)),
    "ap": (Fraction(3), Fraction(5), Fraction(7), Fraction(11)),
    "scaled_main_rows1": (Fraction(4), Fraction(6), Fraction(8), Fraction(10)),
    "mu": (Fraction(2), Fraction(3), Fraction(5), Fraction(7)),
}

MATRIX = {
    "A": (
        (Fraction(2), Fraction(1), Fraction(0), Fraction(1)),
        (Fraction(1), Fraction(0), Fraction(1), Fraction(1)),
    ),
    "D": (
        (Fraction(3), Fraction(2), Fraction(1), Fraction(1)),
        (Fraction(1), Fraction(1), Fraction(0), Fraction(2)),
    ),
    "positive_determinant": (Fraction(1), Fraction(-1), Fraction(2), Fraction(5)),
    "selected_C1": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "selected_C2": (Fraction(-30), Fraction(-135), Fraction(-584), Fraction(-2343)),
    "incoming_rhs1": (Fraction(-30), Fraction(-171), Fraction(-882), Fraction(-4252)),
    "scaled_rhs1": (Fraction(1), Fraction(2), Fraction(3), Fraction(4)),
}


def _matrix_kwargs() -> dict[str, Any]:
    incoming_product = _ordinary_product(HAND_INPUTS["incoming_factor1"], HAND_INPUTS["incoming_m1"])
    scaled_first = _ordinary_neg(
        _ordinary_product(HAND_INPUTS["mu"], HAND_INPUTS["scaled_main_rows1"])
    )
    incoming_first = _ordinary_neg(_ordinary_product(HAND_INPUTS["mu"], incoming_product))
    return {
        "A": MATRIX["A"],
        "D": MATRIX["D"],
        "positive_determinant": MATRIX["positive_determinant"],
        "incoming_rhs": (incoming_first, MATRIX["incoming_rhs1"]),
        "scaled_rhs": (scaled_first, MATRIX["scaled_rhs1"]),
        "ap": HAND_INPUTS["ap"],
        "selected_C1": MATRIX["selected_C1"],
        "selected_C2": MATRIX["selected_C2"],
        "incoming_factor1": HAND_INPUTS["incoming_factor1"],
        "incoming_m1": HAND_INPUTS["incoming_m1"],
        "scaled_main_rows1": HAND_INPUTS["scaled_main_rows1"],
        "mu": HAND_INPUTS["mu"],
    }


def _run_fixture(module: Any, calls: list[dict[str, Any]]) -> dict[str, Any]:
    for name, rows in HAND_INPUTS.items():
        _rows_are_fractional(rows, "hand input " + name)
        if name in {"incoming_factor1", "incoming_m1", "ap", "mu"}:
            _require(all(item != 0 for item in rows), "hand input has a zero history row: " + name)

    incoming_product = _ordinary_product(
        HAND_INPUTS["incoming_factor1"], HAND_INPUTS["incoming_m1"]
    )
    selected_product = _ordinary_product(
        HAND_INPUTS["ap"], HAND_INPUTS["scaled_main_rows1"]
    )
    summed = _ordinary_add(incoming_product, selected_product)
    reduced_expected = _ordinary_neg(_ordinary_product(HAND_INPUTS["mu"], summed))
    _require(incoming_product == (Fraction(2), Fraction(7), Fraction(23), Fraction(72)), "hand incoming product changed")
    _require(selected_product == (Fraction(12), Fraction(38), Fraction(112), Fraction(320)), "hand selected product changed")
    _require(reduced_expected == (Fraction(-28), Fraction(-132), Fraction(-610), Fraction(-2772)), "hand reduced rows changed")

    reduce = _callable(module, "reduce_f_gap_selected_first_row")
    reduction = _invoke(
        calls,
        "reduce_hand_ordinary_Z0_to_Z3",
        reduce,
        HAND_INPUTS["incoming_factor1"],
        HAND_INPUTS["incoming_m1"],
        HAND_INPUTS["ap"],
        HAND_INPUTS["scaled_main_rows1"],
        HAND_INPUTS["mu"],
    )
    _require(_rows(_field(reduction, "incoming_product"), "incoming product") == incoming_product, "incoming product mismatch")
    _require(_rows(_field(reduction, "selected_product"), "selected product") == selected_product, "selected product mismatch")
    _require(_rows(_field(reduction, "reduced"), "reduced") == reduced_expected, "reduced first row mismatch")
    _require(_field(reduction, "coordinate") == "fixed_xi", "coordinate convention changed")
    _assert_uncertified(reduction, "reduction")

    residuals = _callable(module, "matrix_algebra_residuals")
    matrix_args = _matrix_kwargs()
    for name in ("positive_determinant", "selected_C1", "selected_C2", "incoming_rhs1", "scaled_rhs1"):
        _rows_are_fractional(MATRIX[name], "matrix " + name)
        _require(MATRIX[name][3] != 0, "matrix " + name + " omitted its third derivative")
    for name in ("A", "D"):
        rows = MATRIX[name]
        _require(len(rows) == 2, "matrix " + name + " is not a 2-row matrix")
        for index, row in enumerate(rows):
            _rows_are_fractional(row, "matrix " + name + " row " + str(index))
            _require(row[3] != 0, "matrix " + name + " row " + str(index) + " omitted its third derivative")
    for name in ("incoming_rhs", "scaled_rhs"):
        for index, row in enumerate(matrix_args[name]):
            _rows_are_fractional(row, name + " row " + str(index))
            _require(row[3] != 0, name + " row " + str(index) + " omitted its third derivative")

    base = _invoke(calls, "matrix_negative_det_base", residuals, **matrix_args)
    determinant_signed = _ordinary_add(
        _ordinary_product(MATRIX["A"][0], MATRIX["D"][1]),
        _ordinary_neg(_ordinary_product(MATRIX["A"][1], MATRIX["D"][0])),
    )
    _require(determinant_signed == (Fraction(-1), Fraction(1), Fraction(-2), Fraction(-5)), "signed determinant rows changed")
    _require(determinant_signed[0] < 0 and MATRIX["positive_determinant"][0] > 0, "negative determinant sign was not retained")
    _require(_rows(_field(base, "determinant_residual"), "determinant residual") == ZERO, "determinant residual is nonzero")

    incoming_first = _ordinary_neg(
        _ordinary_product(
            HAND_INPUTS["mu"],
            _ordinary_product(HAND_INPUTS["incoming_factor1"], HAND_INPUTS["incoming_m1"]),
        )
    )
    scaled_first = _ordinary_neg(
        _ordinary_product(HAND_INPUTS["mu"], HAND_INPUTS["scaled_main_rows1"])
    )
    matrix_first_target = _ordinary_add(
        incoming_first,
        _ordinary_product(HAND_INPUTS["ap"], scaled_first),
    )
    _require(matrix_first_target == reduced_expected, "RHS first row does not equal raw reduction")
    _require(_rows(_field(base, "incoming_first_rhs_residual"), "incoming RHS residual") == ZERO, "incoming RHS residual is nonzero")
    _require(_rows(_field(base, "scaled_first_rhs_residual"), "scaled RHS residual") == ZERO, "scaled RHS residual is nonzero")
    _require(_rows(_field(base, "matrix_first_row_residual"), "matrix first residual") == ZERO, "matrix first residual is nonzero")
    _require(_rows(_field(base, "matrix_divided_row_residual"), "matrix divided residual") == ZERO, "matrix divided residual is nonzero")
    _require(_rows(_field(base, "raw_first_row_reduction_residual"), "raw reduction residual") == ZERO, "raw reduction residual is nonzero")
    _require(_rows(_field(base, "selected_first_row"), "selected first row") == reduced_expected, "selected first row changed")
    _require(_rows(_field(base, "matrix_first_row_target"), "matrix first target") == reduced_expected, "matrix first target changed")
    _assert_uncertified(base, "matrix base")

    bad_kwargs = dict(matrix_args)
    bad_kwargs["selected_C2"] = (
        Fraction(-29),
        Fraction(-135),
        Fraction(-584),
        Fraction(-2343),
    )
    perturbed = _invoke(calls, "matrix_C2_perturbation_detected", residuals, **bad_kwargs)
    perturbation_residual = _rows(_field(perturbed, "matrix_first_row_residual"), "perturbed matrix residual")
    _require(perturbation_residual == (Fraction(1), Fraction(0), Fraction(1), Fraction(1)), "C2 perturbation was not detected")
    _require(any(item != 0 for item in perturbation_residual), "perturbation residual was silently zero")
    _assert_uncertified(perturbed, "matrix perturbation")

    checks = {
        "ordinary_Z_rows_0_to_3_are_exact": True,
        "incoming_factor1_m1_ap_mu_have_nonzero_rows": True,
        "independent_hand_reduction": True,
        "explicit_third_derivative_formula_used": True,
        "negative_determinant_sign": determinant_signed[0] < 0,
        "both_RHS_rows_checked": True,
        "base_matrix_residuals_zero": True,
        "C2_perturbation_has_nonzero_residual": True,
        "source_identity_false": _field(base, "source_identity_proved") is False,
        "numeric_certification_false": True,
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
        "ordinary_Z_orders": [0, 1, 2, 3],
        "hand_inputs": HAND_INPUTS,
        "hand_products": {
            "incoming_factor1_times_m1": incoming_product,
            "ap_times_scaled_main_rows1": selected_product,
            "sum": summed,
            "reduced": reduced_expected,
        },
        "matrix_example": {
            "A": MATRIX["A"],
            "D": MATRIX["D"],
            "signed_determinant": determinant_signed,
            "positive_determinant": MATRIX["positive_determinant"],
            "incoming_rhs": _matrix_kwargs()["incoming_rhs"],
            "scaled_rhs": _matrix_kwargs()["scaled_rhs"],
            "selected_first_row": reduced_expected,
            "perturbed_C2_matrix_first_residual": perturbation_residual,
            "third_derivative_rows": {
                "A": [row[3] for row in MATRIX["A"]],
                "D": [row[3] for row in MATRIX["D"]],
                "incoming_rhs": [row[3] for row in matrix_args["incoming_rhs"]],
                "scaled_rhs": [row[3] for row in matrix_args["scaled_rhs"]],
                "C1_C2": [MATRIX["selected_C1"][3], MATRIX["selected_C2"][3]],
            },
            "third_derivative_formula": "a3*b0 + 3*a2*b1 + 3*a1*b2 + a0*b3",
        },
        "checks": checks,
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "numeric_certification": False,
        "numeric_certified": False,
        "determinant_positivity_proved": False,
        "error_certified": False,
        "uncertified": True,
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


def _validate_fixture(value: Mapping[str, Any]) -> None:
    _require(value.get("passed") is True and value.get("checks_passed") is True, "synthetic fixture did not pass")
    _require(value.get("synthetic_only") is True, "fixture is not synthetic-only")
    _require(value.get("call_count") == 3, "fixture call count changed")
    calls = value.get("calls")
    _require(isinstance(calls, list) and len(calls) == 3, "fixture call log is incomplete")
    _require(all(item.get("invoked") and item.get("returned") and item.get("passed") for item in calls), "fixture call log has an incomplete call")
    checks = value.get("checks")
    _require(isinstance(checks, Mapping) and all(item is True for item in checks.values()), "fixture contains a failed check")
    _require(not _walk_forbidden_true(value), "fixture made a forbidden source or numeric claim")


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
        "determinant_positivity_proved": False,
        "numeric_certification": False,
        "numeric_certified": False,
        "error_certified": False,
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
        "raw_sha256_authenticated_before_entry": True,
        "leaf_imported": leaf_imported,
        "fixture_invocations": 1 if leaf_imported else 0,
        "leaf_api_invocations": len(calls),
        "source_manifest": _json_safe(source),
        "fixture": _json_safe(fixture),
        "calls": _json_safe(calls),
        "flags": flags,
        **flags,
        "premises": {
            "synthetic_only": True,
            "ordinary_Z0_to_Z3": True,
            "source_identity_proved": False,
            "numeric_certification": False,
            "uncertified": True,
        },
        "error": _json_safe(error),
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
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
            "new_leaf_imported": False,
            "fixture_invocations": 0,
            "leaf_api_invocations": 0,
            "source_identity_proved": False,
            "numeric_certification": False,
            "numeric_certified": False,
            "error_certified": False,
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
                "numeric_certification": False,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
