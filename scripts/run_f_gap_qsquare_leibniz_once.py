"""Run the new F ``q-squared`` Leibniz polynomial fixture once.

This file is the persistence boundary for ``f_gap_qsquare_leibniz_v1``.  It
does not import the historical source graph, an ancestor consumer, or an old
producer.  Before the new leaf is imported it authenticates the raw bytes of
the only two declared inputs (this runner and the new leaf) and records an
exclusive entry receipt.  The leaf is then called for one small, entirely
synthetic polynomial fixture (forward, reverse, and zero-length cases).  The runner independently parses the returned
packet against exact rational expectations; a successful fixture is algebraic
evidence only, not an F-source identity, a certified bound, or numeric-source
admission.

The leaf interface is intentionally small.  It exposes one callable named
``qsquare_leibniz_transport``.  The runner supplies caller-owned scalar
callbacks and a Simpson quadrature callback for each of the three local
orientation cases.  This is one synthetic fixture execution; the transport
call is repeated only across the fixture's forward, reverse, and zero-length
case table.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_qsquare_leibniz_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_gap_qsquare_leibniz_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_gap_qsquare_leibniz.py")
REQUIRED_INPUTS = frozenset({RUNNER.as_posix(), LEAF.as_posix()})
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


# The runner's fixture is data-only metadata for the callback definitions and
# the persisted result.  The expected table below is deliberately separate so
# that the runner does not accept a leaf that merely echoes its own expectation.
SYNTHETIC_FIXTURE: dict[str, Any] = {
    "schema": "f_gap_qsquare_leibniz_polynomial_fixture_v1",
    "synthetic_only": True,
    "coefficient": {"numerator": -1, "denominator": 4},
    "polynomials": {
        "M": {"constant": 1, "z": 1, "xi": 1},
        "Q0": {"z": 1, "xi": 1},
        "Q1": {"constant": 3},
        "lower": {"constant": -1, "z": 1},
        "upper": {"z": 2},
    },
    "cases": {
        "forward": {"z": 2, "lower": 1, "upper": 4, "lower_Z": 1, "upper_Z": 2},
        "reverse": {"z": 2, "lower": 4, "upper": 1, "lower_Z": 2, "upper_Z": 1},
        # l = z - 1 and u = 2 z meet at z = -1.  The integral terms vanish,
        # but the two endpoint derivatives differ, so boundary terms remain.
        "zero_length": {
            "z": -1,
            "lower": -2,
            "upper": -2,
            "lower_Z": 1,
            "upper_Z": 2,
        },
    },
}


# Exact independent expectations.  ``measure_term`` is
# c*integral(M_Z*Q0^2), ``q_gap_term`` is
# -2*c*integral(M*Q0*(Q1 - partial_Z Q0)), and the final two entries are the
# moving-upper and moving-lower boundary terms.  All numbers are Fractions so
# comparison never relies on a floating-point tolerance.
EXPECTED_CASES: dict[str, dict[str, Fraction]] = {
    "forward": {
        "z": Fraction(2),
        "M_Z": Fraction(1),
        "Q_derivative_gap": Fraction(2),
        "lower_endpoint": Fraction(1),
        "upper_endpoint": Fraction(4),
        "lower_endpoint_Z": Fraction(1),
        "upper_endpoint_Z": Fraction(2),
        "integral_M_Z_Q0_squared": Fraction(63),
        "integral_M_Q0_E_Q": Fraction(153),
        "measure_term": Fraction(-63, 4),
        "q_gap_term": Fraction(153, 2),
        "upper_boundary_term": Fraction(-126),
        "lower_boundary_term": Fraction(9),
        "residual": Fraction(-225, 4),
    },
    "reverse": {
        "z": Fraction(2),
        "M_Z": Fraction(1),
        "Q_derivative_gap": Fraction(2),
        "lower_endpoint": Fraction(4),
        "upper_endpoint": Fraction(1),
        "lower_endpoint_Z": Fraction(2),
        "upper_endpoint_Z": Fraction(1),
        "integral_M_Z_Q0_squared": Fraction(-63),
        "integral_M_Q0_E_Q": Fraction(-153),
        "measure_term": Fraction(63, 4),
        "q_gap_term": Fraction(-153, 2),
        "upper_boundary_term": Fraction(-9),
        "lower_boundary_term": Fraction(126),
        "residual": Fraction(225, 4),
    },
    "zero_length": {
        "z": Fraction(-1),
        "M_Z": Fraction(1),
        "Q_derivative_gap": Fraction(2),
        "lower_endpoint": Fraction(-2),
        "upper_endpoint": Fraction(-2),
        "lower_endpoint_Z": Fraction(1),
        "upper_endpoint_Z": Fraction(2),
        "integral_M_Z_Q0_squared": Fraction(0),
        "integral_M_Q0_E_Q": Fraction(0),
        "measure_term": Fraction(0),
        "q_gap_term": Fraction(0),
        "upper_boundary_term": Fraction(9),
        "lower_boundary_term": Fraction(-9, 2),
        "residual": Fraction(9, 2),
    },
}


class RunnerInputError(ValueError):
    """Raised before entry when the frozen one-shot contract is invalid."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    """Hash raw bytes only; the input is not decoded by this helper."""

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
        ).encode("utf-8")
        + b"\n"
    )


def _json_safe(value: Any) -> Any:
    """Make a returned leaf packet persistable without importing helpers."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=repr)
    if isinstance(value, Path):
        return value.as_posix()
    return str(value)


def _write_once(path: Path, value: Mapping[str, Any]) -> str:
    """Write one artifact with O_EXCL semantics and return its raw SHA."""

    raw = _json_bytes(_json_safe(value))
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
    """Authenticate exactly the runner and new leaf, with no extra input."""

    files = source.get("files")
    _require(isinstance(files, list) and len(files) == 2,
             "source manifest must contain exactly two file entries")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_path(item.get("path"), f"source.files[{index}]")
        _require(key not in hashes, "duplicate manifest path: " + key)
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        hashes[key] = actual

    _require(frozenset(hashes) == REQUIRED_INPUTS,
             "manifest must contain exactly runner and new leaf, with no extras")
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    _require(isinstance(raw_path, str) and raw_path,
             "pending brief source checkpoint path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(),
             "pending brief source checkpoint path mismatch")
    expected = checkpoint.get("sha256")
    actual = _sha256_file(SOURCE)
    _require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and actual.lower() == expected.lower(),
        "pending brief source checkpoint SHA mismatch",
    )
    return actual


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    _require(not ENTRY.exists() and not RESULT.exists(),
             "one-entry component already entered")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")

    brief = _load_object(BRIEF, "pending brief")
    source = _load_object(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "brief component changed")
    _require(brief.get("phase") == PHASE, "brief phase changed")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal flags changed",
    )
    source_sha = _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "manifest component changed")
    # These guards are optional metadata for older frozen manifests.  When
    # present they must agree with this source-only, single-entry runner.
    for key, expected in (
        ("one_entry_only", True),
        ("synthetic_only", True),
        ("no_source_graph_traversal", True),
    ):
        if key in source:
            _require(source[key] is expected, "manifest guard changed: " + key)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha


def _import_new_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new q-square Leibniz leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_qsquare_leibniz_v1_new_leaf", path
    )
    _require(spec is not None and spec.loader is not None,
             "new q-square Leibniz leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    # This is the only project module import in the runner, and it occurs only
    # after entered.json has been written.
    spec.loader.exec_module(module)
    return module


def _transport_callable(module: Any) -> Callable[..., Any]:
    function = getattr(module, "qsquare_leibniz_transport", None)
    _require(callable(function), "new leaf qsquare_leibniz_transport is missing")
    return function


def _simpson(function: Callable[[float], float], a: float, b: float) -> float:
    """Oriented Simpson rule; the fixture integrands are quadratic."""

    midpoint = 0.5 * (a + b)
    return (b - a) * (
        function(a) + 4.0 * function(midpoint) + function(b)
    ) / 6.0


def _case_callbacks(
    case: Mapping[str, Any],
) -> tuple[dict[str, Callable[[float], float]], float, float, float, float]:
    z = float(case["z"])

    def M(xi: float) -> float:
        return 1.0 + z + xi

    def M_Z(xi: float) -> float:
        del xi
        return 1.0

    def Q0(xi: float) -> float:
        return z + xi

    def Q1(xi: float) -> float:
        del xi
        return 3.0

    def dZ_Q0(xi: float) -> float:
        del xi
        return 1.0

    return (
        {"M": M, "M_Z": M_Z, "Q0": Q0, "Q1": Q1, "dZ_Q0": dZ_Q0},
        float(case["lower"]),
        float(case["upper"]),
        float(case["lower_Z"]),
        float(case["upper_Z"]),
    )


TRANSPORT_FIELDS = frozenset(
    {
        "term_mz_q0_squared",
        "term_m_q0_q1_minus_dz_q0",
        "term_upper_endpoint",
        "term_lower_endpoint",
        "total",
    }
)


def _invoke_transport(function: Callable[..., Any]) -> dict[str, Any]:
    """Run one synthetic fixture over its three orientation cases."""

    returned: dict[str, Any] = {}
    for case_name, case in SYNTHETIC_FIXTURE["cases"].items():
        callbacks, lower, upper, lower_z, upper_z = _case_callbacks(case)
        terms = function(
            callbacks["M"],
            callbacks["M_Z"],
            callbacks["Q0"],
            callbacks["Q1"],
            callbacks["dZ_Q0"],
            lower,
            upper,
            lower_z,
            upper_z,
            _simpson,
        )
        as_dict = getattr(terms, "as_dict", None)
        _require(callable(as_dict), "transport result has no as_dict()")
        packet = as_dict()
        _require(isinstance(packet, Mapping), "transport as_dict() did not return an object")
        _require(set(packet) == TRANSPORT_FIELDS,
                 "transport result fields changed: expected five signed terms")
        returned[case_name] = dict(packet)
    return returned


def _norm_key(key: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(key).lower())


def _lookup(value: Any, aliases: tuple[str, ...], depth: int = 0) -> Any:
    """Find one field in a small returned case packet, not in source data."""

    if not isinstance(value, Mapping) or depth > 2:
        return None
    wanted = {_norm_key(alias) for alias in aliases}
    for key, item in value.items():
        if _norm_key(key) in wanted:
            return item
    for item in value.values():
        if isinstance(item, Mapping):
            found = _lookup(item, aliases, depth + 1)
            if found is not None:
                return found
    return None


def _exact(value: Any, label: str) -> Fraction:
    """Parse a rational spelling for comparison, without certifying floats.

    A returned ``float`` is converted from its decimal text spelling only so
    that this runner can compare the leaf packet to the hand-derived rational
    table.  This does not recover exact floating-point arithmetic, establish a
    rounding error bound, or turn the result into a certified numerical claim.
    """

    if isinstance(value, bool) or value is None:
        raise RunnerInputError(label + " is not an exact scalar")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float):
        _require(value == value and value not in (float("inf"), float("-inf")),
                 label + " is not finite")
        return Fraction(str(value))
    if isinstance(value, str):
        text = value.strip()
        _require(bool(text), label + " is empty")
        try:
            return Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise RunnerInputError(label + " is not an exact rational: " + text) from exc
    if isinstance(value, Mapping):
        for name in ("exact", "rational", "value", "scalar"):
            if name in value:
                return _exact(value[name], label + "." + name)
        if "numerator" in value and "denominator" in value:
            try:
                return Fraction(int(value["numerator"]), int(value["denominator"]))
            except (TypeError, ValueError, ZeroDivisionError) as exc:
                raise RunnerInputError(label + " has malformed numerator/denominator") from exc
    raise RunnerInputError(label + " is not an exact rational")


TERM_ALIASES: dict[str, tuple[str, ...]] = {
    "measure_term": (
        "term_mz_q0_squared", "measure_term", "M_Z_term",
    ),
    "q_gap_term": (
        "term_m_q0_q1_minus_dz_q0", "q_gap_term", "Q_gap_term",
    ),
    "upper_boundary_term": (
        "term_upper_endpoint", "upper_boundary_term", "upper_term",
    ),
    "lower_boundary_term": (
        "term_lower_endpoint", "lower_boundary_term", "lower_term",
    ),
    "residual": ("total", "residual", "leibniz_residual"),
}


def _validate_transport_results(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "transport calls did not return an object")
    cases = value
    missing = sorted(set(EXPECTED_CASES) - set(cases))
    _require(not missing, "transport results missing cases: " + ", ".join(missing))
    validations: dict[str, Any] = {}
    for case_name, expected in EXPECTED_CASES.items():
        actual_case = cases[case_name]
        checks: dict[str, bool] = {}
        parsed: dict[str, str] = {}
        for field in ("measure_term", "q_gap_term", "upper_boundary_term",
                      "lower_boundary_term", "residual"):
            expected_value = {
                "measure_term": expected["measure_term"],
                "q_gap_term": expected["q_gap_term"],
                "upper_boundary_term": expected["upper_boundary_term"],
                "lower_boundary_term": expected["lower_boundary_term"],
                "residual": expected["residual"],
            }[field]
            raw = _lookup(actual_case, TERM_ALIASES[field])
            _require(raw is not None,
                     f"{case_name} case omitted exact field: {field}")
            actual = _exact(raw, f"{case_name}.{field}")
            checks[field] = actual == expected_value
            _require(checks[field],
                     f"{case_name}.{field} mismatch: expected {expected_value}, got {actual}")
            parsed[field] = str(actual)

        # Independently reconstruct the four Leibniz terms and residual from
        # the parsed leaf values; this catches a self-consistent but wrong
        # residual packet.
        rebuilt = (
            _exact(parsed["measure_term"], "measure_term")
            + _exact(parsed["q_gap_term"], "q_gap_term")
            + _exact(parsed["upper_boundary_term"], "upper_boundary_term")
            + _exact(parsed["lower_boundary_term"], "lower_boundary_term")
        )
        _require(rebuilt == expected["residual"],
                 f"{case_name}.residual does not equal its four parsed terms")
        # These independent values are from the callback polynomial and bounds,
        # not from a leaf-provided expectation packet.
        independent_inputs = {
            key: str(expected[key])
            for key in (
                "z", "M_Z", "Q_derivative_gap", "lower_endpoint",
                "upper_endpoint", "lower_endpoint_Z", "upper_endpoint_Z",
                "integral_M_Z_Q0_squared", "integral_M_Q0_E_Q",
            )
        }
        validations[case_name] = {
            "independent_inputs": independent_inputs,
            "expected_terms": {
                key: str(expected[key])
                for key in (
                    "measure_term", "q_gap_term", "upper_boundary_term",
                    "lower_boundary_term", "residual",
                )
            },
            "parsed": parsed,
            "checks": checks,
            "four_term_sum": str(rebuilt),
        }

    _require(EXPECTED_CASES["forward"]["M_Z"] != 0,
             "forward M_Z unexpectedly vanishes")
    _require(EXPECTED_CASES["forward"]["Q_derivative_gap"] != 0,
             "forward Q derivative gap unexpectedly vanishes")
    _require(
        EXPECTED_CASES["forward"]["lower_endpoint"]
        < EXPECTED_CASES["forward"]["upper_endpoint"],
        "forward endpoints are not ordered",
    )
    _require(
        EXPECTED_CASES["reverse"]["lower_endpoint"]
        > EXPECTED_CASES["reverse"]["upper_endpoint"],
        "reverse endpoints are not reversed",
    )
    _require(
        EXPECTED_CASES["zero_length"]["lower_endpoint"]
        == EXPECTED_CASES["zero_length"]["upper_endpoint"],
        "zero-length endpoints do not coincide",
    )
    _require(
        EXPECTED_CASES["zero_length"]["lower_endpoint_Z"]
        != EXPECTED_CASES["zero_length"]["upper_endpoint_Z"],
        "zero-length endpoint derivatives unexpectedly coincide",
    )
    return {
        "passed": True,
        "comparison_to_hand_derived_rational_table": validations,
        "float_to_fraction_is_decimal_text_comparison_only": True,
        "float_exactness_or_certified_bound": False,
        "nonzero_M_Z": True,
        "nonzero_Q_derivative_gap": True,
        "reverse_bounds_checked": True,
        "zero_length_checked": True,
        "zero_length_endpoint_coincident": True,
        "zero_length_endpoint_derivatives_differ": True,
    }


def _failed_output(
    *,
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    entry_sha: str,
    callable_name: str | None,
    transport_results: Any,
    leaf_import_attempted: bool,
    leaf_imported: bool,
    exc: BaseException,
) -> dict[str, Any]:
    return {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": "failed",
        "source_checkpoint_sha256": source_sha,
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "manifest": _json_safe(source),
        "synthetic_fixture_input": SYNTHETIC_FIXTURE,
        "transport_callable": callable_name,
        "leaf_import_attempted": leaf_import_attempted,
        "leaf_imported": leaf_imported,
        "leaf_import_status": (
            "failed" if leaf_import_attempted and not leaf_imported
            else "succeeded" if leaf_imported else "not_attempted"
        ),
        "fixture_invocations": 1 if callable_name is not None else 0,
        "transport_invocations": len(transport_results) if isinstance(transport_results, Mapping) else 0,
        "source_graph_traversed": False,
        "ancestor_replayed": False,
        "network_called": False,
        "old_code_called": False,
        "source_identity_proved": False,
        "certified": False,
        "numeric_use_admissible": False,
        "runtime_model_metadata": None,
        "transport_results": _json_safe(transport_results),
        "error": {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        },
        "limitations": [
            "This is a synthetic polynomial check only.",
            "Float outputs are parsed from decimal text for rational comparison only; this is not exact floating arithmetic or a certified bound.",
            "No actual F-source identity or certified bound is claimed.",
        ],
    }


def main() -> int:
    # No artifact is consumed and no new leaf is imported until all preflight
    # checks, including both raw file hashes, have passed.
    source, hashes, source_sha = _preflight()
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "manifest_paths_exact": sorted(hashes),
            "new_leaf_imported": False,
            "fixture_invocations": 0,
            "source_graph_traversed": False,
            "ancestor_replayed": False,
            "network_called": False,
            "old_code_called": False,
            "source_identity_proved": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
        },
    )

    callable_name: str | None = None
    transport_results: Any = None
    leaf_import_attempted = False
    leaf_imported = False
    try:
        leaf_import_attempted = True
        module = _import_new_leaf()
        leaf_imported = True
        function = _transport_callable(module)
        callable_name = "qsquare_leibniz_transport"
        transport_results = _invoke_transport(function)
        validation = _validate_transport_results(transport_results)
        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "manifest": _json_safe(source),
            "synthetic_fixture_input": SYNTHETIC_FIXTURE,
            "transport_callable": callable_name,
            "leaf_import_attempted": leaf_import_attempted,
            "leaf_imported": leaf_imported,
            "leaf_import_status": "succeeded",
            "fixture_invocations": 1,
            "transport_invocations": len(transport_results),
            "new_leaf_imported_after_entry": True,
            "source_graph_traversed": False,
            "ancestor_replayed": False,
            "network_called": False,
            "old_code_called": False,
            "source_identity_proved": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
            "transport_results": _json_safe(transport_results),
            "independent_validation": validation,
            "error": None,
            "limitations": [
                "Synthetic polynomial algebra only; no actual F-source identity.",
                "Float outputs are parsed from decimal text for rational comparison only; this is not exact floating arithmetic or a certified bound.",
                "No certified reverse bound, numerical-source admission, or error certificate.",
            ],
        }
    except BaseException as exc:
        output = _failed_output(
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            entry_sha=entry_sha,
            callable_name=callable_name,
            transport_results=transport_results,
            leaf_import_attempted=leaf_import_attempted,
            leaf_imported=leaf_imported,
            exc=exc,
        )

    # Complete success or failure is persisted before the status/SHA is printed.
    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "transport_callable": callable_name,
                "leaf_import_attempted": output.get("leaf_import_attempted", False),
                "leaf_imported": output.get("leaf_imported", False),
                "fixture_invocations": output["fixture_invocations"],
                "transport_invocations": output.get("transport_invocations", 0),
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
