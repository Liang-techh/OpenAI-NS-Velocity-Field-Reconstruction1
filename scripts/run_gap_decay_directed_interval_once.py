"""Run the directed gap-decay interval leaf exactly once.

The runner authenticates only the three frozen source files named by the
component checkpoint, records an exclusive entry before any new leaf call,
and then performs one source-text AST check plus bounded exact-Fraction leaf
calls. The historical producer is never imported or called: its source is
read as text and parsed with :mod:`ast` only.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys
import traceback
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "gap_decay_directed_interval_v1"
PHASE = "gap_decay_directed_interval_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp/gap_decay_directed_interval_v1"
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src/openai_ns_reconstruction/gap_decay_directed_interval.py"
PRODUCER = ROOT / (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py"
)
RUNNER = Path(__file__).resolve()

EXPECTED_FILES = {
    LEAF.resolve(): "src/openai_ns_reconstruction/gap_decay_directed_interval.py",
    RUNNER.resolve(): "scripts/run_gap_decay_directed_interval_once.py",
    PRODUCER.resolve(): (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py"
    ),
}


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _write_new(path: Path, value: Any) -> str:
    _require(not path.exists(), "artifact already exists: " + str(path))
    raw = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _resolve_pinned_path(raw_path: object) -> Path:
    _require(isinstance(raw_path, str) and raw_path, "source file path is missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    return candidate.resolve()


def _verify_source_manifest(frozen: object) -> tuple[str, dict[str, str]]:
    _require(isinstance(frozen, dict), "source checkpoint is not an object")
    _require(
        frozen.get("component") == COMPONENT,
        "source checkpoint component mismatch",
    )
    _require(
        frozen.get("one_entry_only") is True,
        "source checkpoint is not one-entry-only",
    )
    _require(
        frozen.get("schema") == "gap_decay_directed_interval_checkpoint_v1",
        "source checkpoint schema mismatch",
    )
    entries = frozen.get("files")
    _require(
        isinstance(entries, list) and len(entries) == 3,
        "source files must contain exactly three entries",
    )

    seen: set[Path] = set()
    pins: dict[str, str] = {}
    for item in entries:
        _require(isinstance(item, dict), "source file entry is not an object")
        path = _resolve_pinned_path(item.get("path"))
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and len(expected) == 64,
            "source file SHA is invalid",
        )
        _require(
            path in EXPECTED_FILES,
            "source checkpoint contains an unexpected file: " + str(path),
        )
        _require(path not in seen, "source checkpoint contains a duplicate file")
        _require(path.is_file(), "frozen source file is missing: " + str(path))
        actual = _sha(path)
        _require(actual == expected, "frozen source file changed: " + str(path))
        seen.add(path)
        pins[EXPECTED_FILES[path]] = expected

    _require(
        seen == set(EXPECTED_FILES),
        "source checkpoint omits a required raw-SHA file",
    )
    return _sha(SOURCE), pins


def _checkpoint_sha(*briefs: dict[str, Any]) -> str:
    for brief in briefs:
        for key in (f"{COMPONENT}_source_checkpoint", "source_checkpoint"):
            candidate = brief.get(key)
            if not isinstance(candidate, dict):
                continue
            raw_path = candidate.get("path")
            expected = candidate.get("sha256")
            if not isinstance(raw_path, str) or not isinstance(expected, str):
                continue
            if _resolve_pinned_path(raw_path) == SOURCE.resolve():
                return expected
    raise RuntimeError("pending brief source checkpoint for this component is missing")


def _verify_pending_brief(state: object, source_sha256: str) -> dict[str, Any]:
    _require(isinstance(state, dict), "pending brief is not an object")
    _require(state.get("goal_active") is True, "pending brief goal is not active")
    _require(
        state.get("full_goal_complete") is False,
        "pending brief goal is already complete",
    )
    _require(
        state.get("do_not_reenter") is True,
        "pending brief does not forbid re-entry",
    )

    descriptor: dict[str, Any] | None = None
    if state.get("component") == COMPONENT:
        descriptor = state
    else:
        nested = state.get(COMPONENT)
        if isinstance(nested, dict):
            descriptor = nested
        pending = state.get("pending")
        if (
            descriptor is None
            and isinstance(pending, dict)
            and pending.get("component") == COMPONENT
        ):
            descriptor = pending
    _require(descriptor is not None, "pending brief component mismatch")
    _require(descriptor.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        _checkpoint_sha(descriptor, state) == source_sha256,
        "pending brief source checkpoint SHA mismatch",
    )
    return descriptor


def _tree_dump(node: ast.AST) -> str:
    return ast.dump(node, annotate_fields=True, include_attributes=False)


def _ast_binding_check() -> dict[str, Any]:
    """Match the two exact build_graph expression trees without importing it."""

    source_text = PRODUCER.read_text(encoding="utf-8")
    tree = ast.parse(source_text, filename=str(PRODUCER))
    functions = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "build_graph"
    ]
    _require(len(functions) == 1, "producer build_graph definition is not unique")
    build_graph = functions[0]

    expected_p = ast.parse(
        "g.add(g.unary('exp',g.constant(40)),g.constant(11))",
        mode="eval",
    ).body
    expected_mu = ast.parse(
        "g.unary('exp',g.sub(g.unary('log',g.constant('1/1000')),g.mul(g.constant(4),P)))",
        mode="eval",
    ).body

    matches: dict[str, list[ast.Assign]] = {"P": [], "mu": []}
    for statement in build_graph.body:
        if not isinstance(statement, ast.Assign) or len(statement.targets) != 1:
            continue
        target = statement.targets[0]
        if isinstance(target, ast.Name) and target.id in matches:
            matches[target.id].append(statement)

    _require(len(matches["P"]) == 1, "build_graph P assignment is not unique")
    _require(len(matches["mu"]) == 1, "build_graph mu assignment is not unique")
    p_assignment = matches["P"][0]
    mu_assignment = matches["mu"][0]
    _require(
        _tree_dump(p_assignment.value) == _tree_dump(expected_p),
        "build_graph P AST mismatch",
    )
    _require(
        _tree_dump(mu_assignment.value) == _tree_dump(expected_mu),
        "build_graph mu AST mismatch",
    )

    return {
        "status": "matched",
        "producer_path": str(PRODUCER),
        "producer_sha256": _sha(PRODUCER),
        "build_graph_lineno": build_graph.lineno,
        "P": {"lineno": p_assignment.lineno, "exact_tree": _tree_dump(p_assignment.value)},
        "mu": {"lineno": mu_assignment.lineno, "exact_tree": _tree_dump(mu_assignment.value)},
        "producer_imported": False,
        "producer_called": False,
        "graph_called": False,
    }


def _fraction(value: Fraction) -> dict[str, Any]:
    _require(isinstance(value, Fraction), "leaf returned a non-Fraction endpoint")
    return {
        "type": "Fraction",
        "numerator": value.numerator,
        "denominator": value.denominator,
        "text": str(value),
    }


def _record_value(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction(value)
    if isinstance(value, tuple):
        return [_record_value(item) for item in value]
    if isinstance(value, list):
        return [_record_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _record_value(item) for key, item in value.items()}
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise RuntimeError("unsupported non-exact result type: " + type(value).__name__)


def _arguments_record(args: tuple[Any, ...], kwargs: dict[str, Any]) -> dict[str, Any]:
    return {"args": _record_value(args), "kwargs": _record_value(kwargs)}


def _call(
    calls: list[dict[str, Any]],
    label: str,
    function: Callable[..., Any],
    *args: Any,
    **kwargs: Any,
) -> Any:
    record: dict[str, Any] = {
        "label": label,
        "arguments": _arguments_record(args, kwargs),
    }
    try:
        value = function(*args, **kwargs)
    except Exception as exc:
        record.update(
            {
                "status": "unexpected_exception",
                "exception": type(exc).__name__ + ": " + str(exc),
            }
        )
        calls.append(record)
        raise
    record.update({"status": "returned", "output": _record_value(value)})
    calls.append(record)
    return value


def _expected_error(
    calls: list[dict[str, Any]],
    label: str,
    function: Callable[..., Any],
    expected_type: type[BaseException],
    *args: Any,
    **kwargs: Any,
) -> None:
    record: dict[str, Any] = {
        "label": label,
        "arguments": _arguments_record(args, kwargs),
        "expected_exception": expected_type.__name__,
    }
    try:
        value = function(*args, **kwargs)
    except expected_type as exc:
        record.update(
            {
                "status": "expected_exception",
                "exception": type(exc).__name__ + ": " + str(exc),
            }
        )
        calls.append(record)
        return
    except Exception as exc:
        record.update(
            {
                "status": "wrong_exception",
                "exception": type(exc).__name__ + ": " + str(exc),
            }
        )
        calls.append(record)
        raise RuntimeError(label + " raised the wrong exception") from exc
    record.update({"status": "missing_exception", "returned": _record_value(value)})
    calls.append(record)
    raise RuntimeError(label + " did not fail closed")


def _assert_interval(value: Any, *, strict: bool = False) -> tuple[Fraction, Fraction]:
    _require(isinstance(value, tuple) and len(value) == 2, "leaf interval shape mismatch")
    lower, upper = value
    _require(
        isinstance(lower, Fraction) and isinstance(upper, Fraction),
        "leaf interval is not exact",
    )
    _require(
        Fraction(0, 1) <= lower <= upper <= Fraction(1, 1),
        "leaf interval is outside [0, 1]",
    )
    if strict:
        _require(lower < upper, "synthetic interval did not retain a nontrivial enclosure")
    return lower, upper


def _load_leaf() -> Any:
    spec = importlib.util.spec_from_file_location("_gap_decay_directed_interval_once", LEAF)
    _require(spec is not None and spec.loader is not None, "new leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    _require(
        callable(getattr(module, "gap_decay_interval", None)),
        "gap_decay_interval API missing",
    )
    _require(
        callable(getattr(module, "canonical_native_decay_enclosure", None)),
        "canonical API missing",
    )
    return module


def _run_new_leaf(calls: list[dict[str, Any]]) -> dict[str, Any]:
    module = _load_leaf()
    gap_decay_interval = module.gap_decay_interval
    canonical_native_decay_enclosure = module.canonical_native_decay_enclosure

    zero = _call(
        calls,
        "D_zero",
        gap_decay_interval,
        Fraction(0, 1),
        Fraction(0, 1),
        max_terms=8,
        max_bits=16,
        tail_exponent_cap=16,
    )
    _require(zero == (Fraction(0, 1), Fraction(0, 1)), "D=0 result mismatch")

    finite_interval = _call(
        calls,
        "D_fraction_interval",
        gap_decay_interval,
        Fraction(1, 2),
        Fraction(3, 2),
        max_terms=12,
        max_bits=64,
        tail_exponent_cap=64,
    )
    finite_lower, finite_upper = _assert_interval(finite_interval, strict=True)

    large_interval = _call(
        calls,
        "D_large_interval",
        gap_decay_interval,
        Fraction(2000, 1),
        Fraction(2001, 1),
        max_terms=2,
        max_bits=128,
        tail_exponent_cap=128,
    )
    expected_large = (Fraction(1, 1) - Fraction(1, 1 << 128), Fraction(1, 1))
    _require(large_interval == expected_large, "large-D directed bound mismatch")

    canonical_interval = _call(
        calls,
        "canonical_source_bound_bits_128",
        canonical_native_decay_enclosure,
        128,
    )
    _require(canonical_interval == expected_large, "canonical directed bound mismatch")
    _assert_interval(canonical_interval)

    _expected_error(
        calls,
        "budget_max_terms_zero",
        gap_decay_interval,
        ValueError,
        Fraction(1, 1),
        Fraction(2, 1),
        max_terms=0,
        max_bits=16,
        tail_exponent_cap=16,
    )
    _expected_error(
        calls,
        "budget_max_bits_above_hard_cap",
        gap_decay_interval,
        ValueError,
        Fraction(1, 1),
        Fraction(2, 1),
        max_terms=8,
        max_bits=4097,
        tail_exponent_cap=128,
    )
    _expected_error(
        calls,
        "budget_tail_cap_above_hard_cap",
        gap_decay_interval,
        ValueError,
        Fraction(1, 1),
        Fraction(2, 1),
        max_terms=8,
        max_bits=128,
        tail_exponent_cap=4097,
    )
    _expected_error(
        calls,
        "reversed_D_interval",
        gap_decay_interval,
        ValueError,
        Fraction(2, 1),
        Fraction(1, 1),
        max_terms=8,
        max_bits=16,
        tail_exponent_cap=16,
    )

    return {
        "canonical_source_premise": {
            "mu_formula": "mu=(1/1000)*exp(-4*(exp(40)+11))",
            "mu_positive": True,
            "mu_less_than_one_over_1000": True,
            "D_formula": "D=2/mu-4",
            "D_strict_lower_bound": "D>1996",
        },
        "canonical_gap_interval": {
            "kind": "directed_closed_interval",
            "not_a_function_value": True,
            "lower": _fraction(canonical_interval[0]),
            "upper": _fraction(canonical_interval[1]),
        },
        "synthetic": {
            "D_zero_exact": _record_value(zero),
            "D_fraction_interval": {
                "input": {
                    "D_lo": _fraction(Fraction(1, 2)),
                    "D_hi": _fraction(Fraction(3, 2)),
                },
                "output": _record_value(finite_interval),
                "strict_enclosure": finite_lower < finite_upper,
            },
            "D_large_interval_exact": _record_value(large_interval),
            "budget_and_order_rejections_recorded": True,
        },
        "F_ownership": False,
        "midpoint_used": False,
        "floating_point_used": False,
        "source_values_selected": False,
        "ancestor_math_calls": 0,
        "graph_calls": 0,
    }


def main() -> int:
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "frozen source checkpoint is missing")
    _require(LEAF.is_file(), "new gap-decay leaf is missing")
    _require(PRODUCER.is_file(), "bound producer source is missing")
    _require(
        not ENTRY.exists() and not RESULT.exists(),
        "component already entered; never reenter",
    )

    frozen = json.loads(SOURCE.read_bytes())
    source_sha256, source_pins = _verify_source_manifest(frozen)
    brief = json.loads(BRIEF.read_bytes())
    pending = _verify_pending_brief(brief, source_sha256)

    entry_sha256 = _write_new(
        ENTRY,
        {
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_sha256": source_sha256,
            "source_pins": source_pins,
            "pending_brief_phase": pending.get("phase"),
            "ast_calls_started": 0,
            "new_leaf_calls_started": 0,
            "ancestor_math_calls": 0,
            "producer_imported": False,
            "producer_called": False,
            "graph_calls": 0,
            "runtime_model_metadata": None,
        },
    )

    calls: list[dict[str, Any]] = []
    try:
        ast_result = _ast_binding_check()
        leaf_result = _run_new_leaf(calls)
        output: dict[str, Any] = {
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "source_pins": source_pins,
            "ast_calls": 1,
            "ast_result": ast_result,
            "new_leaf_calls": len(calls),
            "consumer_calls": 1,
            "ancestor_math_calls": 0,
            "F_ownership": False,
            "source_values_selected": False,
            "producer_imported": False,
            "producer_called": False,
            "graph_calls": 0,
            "runtime_model_metadata": None,
            "call_log": calls,
            "result": leaf_result,
        }
    except Exception as exc:
        output = {
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "source_pins": source_pins,
            "ast_calls": 1 if "ast_result" in locals() else 0,
            "new_leaf_calls": len(calls),
            "consumer_calls": 1,
            "ancestor_math_calls": 0,
            "F_ownership": False,
            "source_values_selected": False,
            "producer_imported": False,
            "producer_called": False,
            "graph_calls": 0,
            "runtime_model_metadata": None,
            "call_log": calls,
            "exception": type(exc).__name__ + ": " + str(exc),
            "traceback": traceback.format_exc(),
        }

    result_sha256 = _write_new(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha256,
                "consumer_calls": 1,
                "new_leaf_calls": len(calls),
                "ancestor_math_calls": 0,
                "runtime_model_metadata": None,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
