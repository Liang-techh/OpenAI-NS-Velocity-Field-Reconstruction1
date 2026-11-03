"""Run the reserved same-graph C1 matrix leaf exactly once.

The frozen ``tmp/reserved_matrix_c1_source_v2/source.json`` is the only input
manifest.  Every declared repository-relative input is authenticated by its
raw SHA-256 before the new matrix module is imported.  The checks below use a
single mock graph whose handles carry exact :class:`fractions.Fraction`
values; no archive, ancestor graph, fixture, or old producer is imported.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "reserved_matrix_c1_source_v2"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = Path("src/openai_ns_reconstruction/reserved_matrix_c1_action_v1.py")
RUNNER = Path("scripts/run_reserved_matrix_c1_source_v2_once.py")
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry source contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
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


def _load_json(raw: bytes, label: str) -> dict[str, Any]:
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


def _write_exclusive(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
    try:
        with path.open("xb") as stream:
            stream.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry receipt already exists: " + str(path)) from exc
    return _sha256_bytes(raw)


def _repository_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(),
             label + " is missing or escapes repository")
    return key, path


def _manifest_records(source: Mapping[str, Any]) -> list[tuple[Any, Any, str]]:
    records: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            expected = item.get("sha256", item.get("raw_sha256"))
            records.append((item.get("path"), expected, f"source.files[{index}]"))
    for field in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            records.extend((key, value, f"source.{field}") for key, value in table.items())
    _require(records, "source manifest contains no raw SHA entries")
    return records


def _authenticate_inputs(source: Mapping[str, Any]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for raw_path, expected, label in _manifest_records(source):
        key, path = _repository_path(raw_path, label)
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        expected = expected.lower()
        prior = hashes.get(key)
        if prior is not None:
            _require(prior == expected, "conflicting SHA entries for " + key)
            continue
        actual = _sha256_file(path)
        _require(actual.lower() == expected, "raw SHA changed: " + key)
        hashes[key] = actual
    _require(LEAF.as_posix() in hashes, "new matrix module is not pinned")
    _require(RUNNER.as_posix() in hashes, "new runner is not pinned")
    return hashes


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(SOURCE.is_file(), "source manifest is missing")
    source_raw = SOURCE.read_bytes()
    source = _load_json(source_raw, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component changed")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")
    hashes = _authenticate_inputs(source)
    return source, hashes, _sha256_bytes(source_raw)


class MockGraph:
    """One live owner plus exact Fraction-valued source handles."""

    def __init__(self) -> None:
        self._next_node = 1

    def handle(self, value: Fraction | int, *, node: int | None = None) -> "MockHandle":
        if node is None:
            node = self._next_node
            self._next_node += 1
        else:
            self._next_node = max(self._next_node, node + 1)
        return MockHandle(self, Fraction(value), node)

    def _check(self, handle: "MockHandle") -> None:
        if not isinstance(handle, MockHandle) or handle.graph is not self:
            raise ValueError("mock operation received a foreign graph handle")

    def add(self, left: "MockHandle", right: "MockHandle") -> "MockHandle":
        self._check(left)
        self._check(right)
        return self.handle(left.value + right.value)

    def mul(self, left: "MockHandle", right: "MockHandle") -> "MockHandle":
        self._check(left)
        self._check(right)
        return self.handle(left.value * right.value)

    def neg(self, value: "MockHandle") -> "MockHandle":
        self._check(value)
        return self.handle(-value.value)


class MockHandle:
    def __init__(self, graph: MockGraph, value: Fraction, node: int) -> None:
        self.graph = graph
        self.value = value
        self.node = node


class MockPair:
    def __init__(self, value: MockHandle, derivative: MockHandle) -> None:
        self.value = value
        self.Z = derivative


class ForgedScalarGraphHandle:
    """Looks named but has no non-scalar live graph owner."""

    def __init__(self, value: Fraction, node: int) -> None:
        self.graph = "graph-a"
        self.value = value
        self.node = node


def _fraction(value: Any) -> Fraction:
    if not isinstance(value, Fraction):
        raise RunnerInputError("matrix result is not an exact Fraction")
    return value


def _values(handles: Any) -> tuple[Fraction, ...]:
    return tuple(_fraction(handle.value) for handle in handles)


def _fraction_json(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _vector_json(values: Any) -> list[dict[str, int]]:
    return [_fraction_json(_fraction(value)) for value in values]


def _pair_json(pair: Any) -> dict[str, dict[str, int]]:
    return {
        "value": _fraction_json(_fraction(pair.value.value)),
        "Z": _fraction_json(_fraction(pair.Z.value)),
    }


def _record(
    checks: dict[str, Any],
    name: str,
    passed: bool,
    *,
    expected: Any = None,
    actual: Any = None,
    error: str | None = None,
) -> None:
    item: dict[str, Any] = {"passed": bool(passed)}
    if expected is not None:
        item["expected"] = expected
    if actual is not None:
        item["actual"] = actual
    if error is not None:
        item["error"] = error
    checks[name] = item
    if not passed:
        raise RunnerInputError("partial check failed: " + name)


def _build_inputs(graph: MockGraph) -> tuple[Any, Any, Any, Any, Any, Any]:
    matrix_values = (
        (1, 1, 0),
        (0, 1, 1),
        (0, 0, 1),
    )
    matrix_derivatives = (
        (Fraction(0), Fraction(1, 2), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1, 3)),
        (Fraction(1, 4), Fraction(0), Fraction(0)),
    )
    inverse_values = (
        (1, -1, 1),
        (0, 1, -1),
        (0, 0, 1),
    )
    y_values = (1, 2, 3)
    y_derivatives = (Fraction(1, 2), Fraction(-1), Fraction(2))
    matrix = tuple(tuple(graph.handle(value) for value in row) for row in matrix_values)
    matrix_Z = tuple(tuple(graph.handle(value) for value in row) for row in matrix_derivatives)
    inverse = tuple(tuple(graph.handle(value) for value in row) for row in inverse_values)
    y = tuple(graph.handle(value) for value in y_values)
    y_Z = tuple(graph.handle(value) for value in y_derivatives)
    return matrix, matrix_Z, inverse, y, y_Z, matrix_derivatives


def _run_matrix_checks(module: Any, checks: dict[str, Any]) -> dict[str, Any]:
    graph = MockGraph()
    matrix, matrix_Z, inverse, y, y_Z, _ = _build_inputs(graph)
    operations = dict(add=graph.add, mul=graph.mul, neg=graph.neg)

    action = module.inverse_c1_action(
        matrix,
        inverse,
        y,
        y_Z,
        matrix_Z=matrix_Z,
        pair_factory=MockPair,
        **operations,
    )
    beta = _values(action.beta)
    beta_Z = _values(action.beta_Z)
    _record(
        checks,
        "beta",
        beta == (Fraction(2), Fraction(-1), Fraction(3)),
        expected=_vector_json((Fraction(2), Fraction(-1), Fraction(3))),
        actual=_vector_json(beta),
    )
    _record(
        checks,
        "beta_Z",
        beta_Z == (Fraction(9, 2), Fraction(-7, 2), Fraction(3, 2)),
        expected=_vector_json((Fraction(9, 2), Fraction(-7, 2), Fraction(3, 2))),
        actual=_vector_json(beta_Z),
    )
    _record(
        checks,
        "same_graph_report",
        action.report.same_graph_handles
        and action.report.source_handles_bound
        and action.report.caller_inverse_assumed
        and not action.report.inverse_identity_verified,
        expected=True,
        actual={
            "same_graph_handles": action.report.same_graph_handles,
            "source_handles_bound": action.report.source_handles_bound,
            "caller_inverse_assumed": action.report.caller_inverse_assumed,
            "inverse_identity_verified": action.report.inverse_identity_verified,
        },
    )

    full = module.full_third_row_c1_action(
        matrix,
        inverse,
        y,
        y_Z,
        matrix_Z=matrix_Z,
        complete_m15_row=matrix[2],
        complete_m15_row_Z=matrix_Z[2],
        pair_factory=MockPair,
        require_full_row_identity=True,
        **operations,
    )
    full_value = _fraction(full.response.value.value)
    full_Z = _fraction(full.response.Z.value)
    _record(
        checks,
        "full_row",
        full_value == Fraction(3),
        expected=_fraction_json(Fraction(3)),
        actual=_fraction_json(full_value),
    )
    _record(
        checks,
        "full_row_Z",
        full_Z == Fraction(2),
        expected=_fraction_json(Fraction(2)),
        actual=_fraction_json(full_Z),
    )
    _record(
        checks,
        "full_row_identity_conditionally_available",
        full.report.third_row_c1_identity_checked
        and full.report.full_row_identity_conditionally_available
        and full.report.conditional_response_row == (0, 0, 1)
        and full.report.conditional_response_row_Z == (0, 0, 0),
        expected=True,
        actual={
            "third_row_c1_identity_checked": full.report.third_row_c1_identity_checked,
            "full_row_identity_conditionally_available": full.report.full_row_identity_conditionally_available,
            "conditional_response_row": full.report.conditional_response_row,
            "conditional_response_row_Z": full.report.conditional_response_row_Z,
        },
    )

    zero_matrix_Z = tuple(
        tuple(graph.handle(Fraction(0)) for _ in range(3)) for _ in range(3)
    )
    matrix_pairs = tuple(
        tuple(MockPair(matrix[i][j], zero_matrix_Z[i][j]) for j in range(3))
        for i in range(3)
    )
    y_pairs = tuple(MockPair(y[i], y_Z[i]) for i in range(3))
    omitted = module.inverse_c1_action(
        matrix_pairs,
        inverse,
        y_pairs,
        pair_factory=MockPair,
        **operations,
    )
    omitted_beta_Z = _values(omitted.beta_Z)
    expected_omitted = (Fraction(7, 2), Fraction(-3), Fraction(2))
    _record(
        checks,
        "omitted_matrix_Z_is_different",
        omitted_beta_Z == expected_omitted and omitted_beta_Z != beta_Z,
        expected={
            "omitted_beta_Z": _vector_json(expected_omitted),
            "different_from": _vector_json(beta_Z),
        },
        actual=_vector_json(omitted_beta_Z),
    )

    other_graph = MockGraph()
    mixed_matrix = [list(row) for row in matrix]
    mixed_matrix[0][0] = other_graph.handle(1)
    try:
        module.inverse_c1_action(
            mixed_matrix,
            inverse,
            y,
            y_Z,
            matrix_Z=matrix_Z,
            pair_factory=MockPair,
            **operations,
        )
    except Exception as exc:
        mixed_error = type(exc).__name__
    else:
        mixed_error = None
    _record(
        checks,
        "reject_mixed_graph_owner",
        mixed_error is not None,
        expected="exception",
        actual=mixed_error,
    )

    forged_matrix = [list(row) for row in matrix]
    forged_matrix[0][0] = ForgedScalarGraphHandle(Fraction(1), 1)
    try:
        module.inverse_c1_action(
            forged_matrix,
            inverse,
            y,
            y_Z,
            matrix_Z=matrix_Z,
            pair_factory=MockPair,
            **operations,
        )
    except Exception as exc:
        forged_error = type(exc).__name__
    else:
        forged_error = None
    _record(
        checks,
        "reject_forged_scalar_graph_id",
        forged_error is not None,
        expected="exception",
        actual=forged_error,
    )

    wrong_third = list(matrix[2])
    wrong_third[2] = graph.handle(1, node=10_000_003)
    try:
        module.full_third_row_c1_action(
            matrix,
            inverse,
            y,
            y_Z,
            matrix_Z=matrix_Z,
            complete_m15_row=tuple(wrong_third),
            complete_m15_row_Z=matrix_Z[2],
            pair_factory=MockPair,
            require_full_row_identity=True,
            **operations,
        )
    except Exception as exc:
        third_row_error = type(exc).__name__
    else:
        third_row_error = None
    _record(
        checks,
        "reject_different_third_row_node",
        third_row_error is not None,
        expected="exception",
        actual=third_row_error,
    )

    return {
        "beta": _vector_json(beta),
        "beta_Z": _vector_json(beta_Z),
        "full_row": _pair_json(full.response),
        "omitted_matrix_Z_beta_Z": _vector_json(omitted_beta_Z),
        "negative_case_exception_types": {
            "mixed_graph_owner": mixed_error,
            "forged_scalar_graph_id": forged_error,
            "different_third_row_node": third_row_error,
        },
    }


def _import_leaf(hashes: Mapping[str, str]) -> Any:
    path = (ROOT / LEAF).resolve()
    _require(_sha256_file(path).lower() == hashes[LEAF.as_posix()].lower(),
             "new matrix module changed after entry")
    module_name = "reserved_matrix_c1_action_v1_new_leaf"
    _require(module_name not in sys.modules, "new matrix module already imported")
    spec = importlib.util.spec_from_file_location(module_name, path)
    _require(spec is not None and spec.loader is not None, "new matrix module load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    source, hashes, source_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "component": COMPONENT,
            "entered": True,
            "source_sha256": source_sha,
            "input_count": len(hashes),
            "input_sha256": dict(sorted(hashes.items())),
        },
    )
    partial_checks: dict[str, Any] = {}
    try:
        module = _import_leaf(hashes)
        _require(getattr(module, "SCHEMA", None) == "reserved_matrix_c1_action_v1",
                 "new matrix module schema changed")
        outputs = _run_matrix_checks(module, partial_checks)
        result: dict[str, Any] = {
            "component": COMPONENT,
            "status": "returned",
            "entry_sha256": entry_sha,
            "source_sha256": source_sha,
            "raw_inputs_authenticated": True,
            "old_math_or_fixture_imported": False,
            "ancestor_graph_read": False,
            "module_imported_after_entry": True,
            "partial_checks": partial_checks,
            "outputs": outputs,
            "traceback": None,
        }
    except BaseException:
        result = {
            "component": COMPONENT,
            "status": "failed",
            "entry_sha256": entry_sha,
            "source_sha256": source_sha,
            "raw_inputs_authenticated": True,
            "old_math_or_fixture_imported": False,
            "ancestor_graph_read": False,
            "module_imported_after_entry": True,
            "partial_checks": partial_checks,
            "traceback": traceback.format_exc(),
        }
    result_sha = _write_exclusive(RESULT, result)
    print(json.dumps({"status": result["status"], "result_sha256": result_sha}, sort_keys=True))


if __name__ == "__main__":
    main()
