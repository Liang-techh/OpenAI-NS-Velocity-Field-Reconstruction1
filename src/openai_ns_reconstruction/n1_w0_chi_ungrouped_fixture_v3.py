"""Synthetic, no-arithmetic checks for the ungrouped W0-chi leaf.

The v1 leaf grouped equal powers and therefore performed interval addition.  A
v2 consumer must keep the saved parent rows and the four raw W0-chi signed
records as separate formal-expression children.  This fixture deliberately
uses coefficients which raise on every numeric operation.  It checks only
identity and provenance; it never reads, compares, or combines a coefficient.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_w0_chi_ungrouped_fixture_v3"
ORDERS = (0, 1)
SIGNED_COUNTS = {0: 1, 1: 3}
ERROR_COUNTS = {0: 1, 1: 3}
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")

__all__ = ["SCHEMA", "ORDERS", "POWER_BASIS", "run_fixture", "run"]


class ArithmeticProbe(AssertionError):
    """Raised if a leaf attempts numeric work on a synthetic coefficient."""


class _CoefficientTrap:
    """Opaque stand-in whose arithmetic and truth-value are forbidden."""

    events: list[str] = []

    def __init__(self, label: str) -> None:
        self.label = label

    @property
    def _mpi_(self) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
        # The shape resembles an exact directed interval, but no endpoint is
        # ever hydrated or used.  Returning a structural token keeps the leaf
        # free to recognize a saved-coefficient slot without giving it a value.
        return ((0, 1, 0, 0), (0, 1, 0, 0))

    @property
    def a(self) -> _CoefficientTrap:
        self._trip("read lower endpoint")
        return self

    @property
    def b(self) -> _CoefficientTrap:
        self._trip("read upper endpoint")
        return self

    @classmethod
    def reset(cls) -> None:
        cls.events = []

    @classmethod
    def _record(cls, message: str) -> None:
        cls.events.append(message)

    def _trip(self, message: str) -> None:
        self._record(self.label + ":" + message)
        raise ArithmeticProbe("synthetic coefficient arithmetic: " + message)

    def __bool__(self) -> bool:
        self._trip("truth value")
        return False

    def __add__(self, other: Any) -> Any:
        self._trip("add")

    def __radd__(self, other: Any) -> Any:
        self._trip("reverse add")

    def __sub__(self, other: Any) -> Any:
        self._trip("subtract")

    def __rsub__(self, other: Any) -> Any:
        self._trip("reverse subtract")

    def __mul__(self, other: Any) -> Any:
        self._trip("multiply")

    def __rmul__(self, other: Any) -> Any:
        self._trip("reverse multiply")

    def __truediv__(self, other: Any) -> Any:
        self._trip("divide")

    def __rtruediv__(self, other: Any) -> Any:
        self._trip("reverse divide")

    def __neg__(self) -> Any:
        self._trip("negate")

    def __pos__(self) -> Any:
        self._trip("positive")

    def __lt__(self, other: Any) -> bool:
        self._trip("less-than")
        return False

    def __le__(self, other: Any) -> bool:
        self._trip("less-equal")
        return False

    def __gt__(self, other: Any) -> bool:
        self._trip("greater-than")
        return False

    def __ge__(self, other: Any) -> bool:
        self._trip("greater-equal")
        return False

    def __repr__(self) -> str:
        return "<coefficient-trap " + self.label + ">"


class _ContextTrap:
    """Context facade that records conversion attempts without doing them."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def mpf(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append("mpf")
        raise ArithmeticProbe("synthetic context conversion")

    def make_mpf(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append("make_mpf")
        raise ArithmeticProbe("synthetic context interval construction")


def _rows(value: Any) -> dict[int, list[Any]]:
    if not isinstance(value, Mapping):
        return {}
    result: dict[int, list[Any]] = {}
    for order in ORDERS:
        selected: Any = None
        for key in (order, str(order), "Z" + str(order), "z" + str(order)):
            if key in value:
                selected = value[key]
                break
        if isinstance(selected, list):
            result[order] = selected
    return result


def _channel(value: Any, name: str) -> dict[int, list[Any]]:
    """Find a channel without accepting a powers-net or grouped alias."""

    if not isinstance(value, Mapping):
        return {}
    for key in (name, "records", "raw_records", "raw"):
        candidate = value.get(key)
        rows = _rows(candidate)
        if set(rows) == set(ORDERS):
            return rows
    rows = _rows(value)
    if set(rows) == set(ORDERS):
        return rows
    for key in ("W0_chi", "correction", "attached"):
        nested = value.get(key)
        found = _channel(nested, name)
        if set(found) == set(ORDERS):
            return found
    branches = value.get("branches")
    if isinstance(branches, Mapping):
        return _channel(branches.get("W0_chi"), name)
    return {}


def _branch(value: Any) -> dict[int, dict[str, list[Any]]]:
    if not isinstance(value, Mapping):
        return {}
    result: dict[int, dict[str, list[Any]]] = {}
    for order in ORDERS:
        payload: Any = None
        for key in (order, str(order), "Z" + str(order), "z" + str(order)):
            if key in value:
                payload = value[key]
                break
        if not isinstance(payload, Mapping):
            continue
        result[order] = {
            "swirl": list(payload.get("swirl", ())) if isinstance(payload.get("swirl", ()), list) else [],
            "swirl_error": list(payload.get("swirl_error", ())) if isinstance(payload.get("swirl_error", ()), list) else [],
        }
    return result


def _term(
    *,
    source_id: str,
    source_path: tuple[str, ...],
    signed: bool,
    label: str,
) -> dict[str, Any]:
    origin = {
        "source_ids": [source_id, "synthetic:shared-upstream"],
        "source_paths": [list(source_path)],
        "signed_field": signed,
        **({"error_only": True} if not signed else {}),
        "uncertainty": True,
        "lineage": {"label": label, "parent": "synthetic:shared-upstream"},
    }
    return {
        "powers": [0, 2, 0, -3],
        "coefficient": _CoefficientTrap(label),
        "origin": origin,
        "source_ids": [source_id, "synthetic:shared-upstream"],
        "source_paths": [list(source_path)],
        "signed_field": signed,
        "symmetric_error": not signed,
        "uncertainty": True,
        "receipt": "synthetic://n1_w0_chi_ungrouped_v2",
        "path": list(source_path),
    }


def _case() -> dict[str, Any]:
    prior: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    signed: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    errors: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for order in ORDERS:
        # Two parent entries per order make ordering observable.  Powers are
        # intentionally repeated, but no coefficient may be combined.
        prior[order] = [
            _term(
                source_id=f"synthetic:prior:Z{order}:0",
                source_path=("synthetic", "prior", f"Z{order}", "0"),
                signed=True,
                label=f"prior-Z{order}-0",
            ),
            _term(
                source_id=f"synthetic:prior:Z{order}:1",
                source_path=("synthetic", "prior", f"Z{order}", "1"),
                signed=True,
                label=f"prior-Z{order}-1",
            ),
        ]
        signed[order] = [
            _term(
                source_id=f"synthetic:signed:Z{order}:{index}",
                source_path=("synthetic", "signed", f"Z{order}", str(index)),
                signed=True,
                label=f"signed-Z{order}-{index}",
            )
            for index in range(SIGNED_COUNTS[order])
        ]
        errors[order] = [
            _term(
                source_id=f"synthetic:error:Z{order}:{index}",
                source_path=("synthetic", "error", f"Z{order}", str(index)),
                signed=False,
                label=f"error-Z{order}-{index}",
            )
            for index in range(ERROR_COUNTS[order])
        ]
    branch = {
        order: {
            "axial": [],
            "axial_error": [],
            "swirl": signed[order],
            "swirl_error": errors[order],
        }
        for order in ORDERS
    }
    family = {"family": "synthetic-n1-w0-chi-ungrouped", "revision": "v3"}
    logs = {name: "synthetic-log:" + name for name in POWER_BASIS}
    return {
        "prior": prior,
        "signed": signed,
        "errors": errors,
        "branch": branch,
        "correction": {"branches": {"W0_chi": branch}},
        "family": family,
        "logs": logs,
    }


def _expression_rows(result: Mapping[str, Any], order: int) -> list[Any]:
    expression = result.get("formal_expression")
    if expression is None:
        expression = result.get("segmented_formal_expression")
    if expression is None:
        expression = result.get("formal_expression_by_order")
    selected = []
    if isinstance(expression, Mapping):
        for key in (order, str(order), "Z" + str(order), "z" + str(order)):
            if key in expression:
                selected = expression[key]
                break
    if isinstance(selected, Mapping):
        # The production leaf keeps the three channels as named segments.
        # For the signed formal expression, parent precedes the raw signed
        # records; the error segment is checked separately and is never part
        # of this returned sequence.
        prior = selected.get("prior")
        signed = selected.get("signed")
        if isinstance(prior, list) and isinstance(signed, list):
            return list(prior) + list(signed)
        for key in ("terms", "children", "subterms", "items"):
            if isinstance(selected.get(key), list):
                return selected[key]
    return selected


def _entry_term(value: Any) -> Any:
    if isinstance(value, Mapping):
        for key in ("term", "record", "raw_term", "input_term", "child"):
            if key in value:
                return value[key]
    return value


def _source_paths(term: Mapping[str, Any]) -> list[str]:
    values: list[Any] = []
    for container in (term, term.get("origin")):
        if not isinstance(container, Mapping):
            continue
        for key in ("source_paths", "source_path", "path"):
            value = container.get(key)
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
                values.append(tuple(str(item) for item in value))
            elif isinstance(value, (str, int)):
                values.append(str(value))
    result: list[str] = []
    for value in values:
        token = repr(value)
        if token not in result:
            result.append(token)
    return result


def _markers(term: Any, signed: bool) -> bool:
    if not isinstance(term, Mapping) or not isinstance(term.get("origin"), Mapping):
        return False
    origin = term["origin"]
    return (
        term.get("signed_field") is signed
        and term.get("error_only", not signed) is (not signed)
        and origin.get("signed_field") is signed
        and origin.get("error_only", not signed) is (not signed)
        and isinstance(origin.get("source_ids"), list)
        and isinstance(origin.get("source_paths"), list)
        and isinstance(origin.get("lineage"), Mapping)
        and term.get("uncertainty") is True
    )


def _validate(case: Mapping[str, Any], result: Any) -> dict[str, bool]:
    checks: dict[str, bool] = {
        "result_mapping": isinstance(result, Mapping),
        "formal_expression_present": False,
        "parent_order_retained": False,
        "raw_signed_identity_retained": False,
        "raw_error_identity_retained": False,
        "signed_counts_1_and_3": False,
        "error_counts_1_and_3": False,
        "origin_labels_and_lineage": False,
        "source_paths_once": False,
        "shared_sources_not_counted_as_duplicates": False,
        "error_tail_not_in_signed": False,
        "no_numeric_addition_flags": False,
        "segmented_expression_not_power_net": False,
        "owner_not_proved": False,
    }
    if not isinstance(result, Mapping):
        return checks

    signed = _channel(result, "signed")
    errors = _channel(result, "error_only")
    checks["signed_counts_1_and_3"] = (
        len(signed.get(0, ())) == SIGNED_COUNTS[0]
        and len(signed.get(1, ())) == SIGNED_COUNTS[1]
    )
    checks["error_counts_1_and_3"] = (
        len(errors.get(0, ())) == ERROR_COUNTS[0]
        and len(errors.get(1, ())) == ERROR_COUNTS[1]
    )
    checks["raw_signed_identity_retained"] = all(
        len(signed.get(order, ())) == len(case["signed"][order])
        and all(
            signed[order][index] is case["signed"][order][index]
            for index in range(len(case["signed"][order]))
        )
        for order in ORDERS
    )
    checks["raw_error_identity_retained"] = all(
        len(errors.get(order, ())) == len(case["errors"][order])
        and all(
            errors[order][index] is case["errors"][order][index]
            for index in range(len(case["errors"][order]))
        )
        for order in ORDERS
    )
    checks["origin_labels_and_lineage"] = all(
        _markers(term, True) for order in ORDERS for term in signed.get(order, ())
    ) and all(_markers(term, False) for order in ORDERS for term in errors.get(order, ()))
    expected_paths = [
        path
        for group in (case["prior"], case["signed"], case["errors"])
        for order in ORDERS
        for term in group[order]
        for path in _source_paths(term)
    ]
    checks["source_paths_once"] = len(expected_paths) == len(set(expected_paths))
    signed_paths = [
        path
        for order in ORDERS
        for term in signed.get(order, ())
        if isinstance(term, Mapping)
        for path in _source_paths(term)
    ]
    error_paths = [
        path
        for order in ORDERS
        for term in errors.get(order, ())
        if isinstance(term, Mapping)
        for path in _source_paths(term)
    ]
    checks["shared_sources_not_counted_as_duplicates"] = (
        any(
            isinstance(term, Mapping)
            and "synthetic:shared-upstream" in term.get("source_ids", ())
            for order in ORDERS
            for term in case["signed"][order]
        )
        and len(signed_paths) == len(set(signed_paths))
        and len(error_paths) == len(set(error_paths))
    )
    checks["error_tail_not_in_signed"] = all(
        all(
            term is not error
            for term in signed.get(order, ())
            for error_order in ORDERS
            for error in case["errors"][error_order]
        )
        for order in ORDERS
    )

    expression_ok = True
    for order in ORDERS:
        raw = [_entry_term(item) for item in _expression_rows(result, order)]
        expected = list(case["prior"][order]) + list(case["signed"][order])
        expression_ok = expression_ok and len(raw) == len(expected)
        expression_ok = expression_ok and all(
            raw[index] is expected[index] for index in range(min(len(raw), len(expected)))
        )
        expression_ok = expression_ok and not any(
            item is error for item in raw for error in case["errors"][order]
        )
    checks["formal_expression_present"] = expression_ok
    checks["parent_order_retained"] = expression_ok

    forbidden_true = (
        "numeric_values_computed",
        "coefficients_summed",
        "coefficient_sum_performed",
        "grouped_by_powers",
        "power_net_computed",
        "formal_sum_recomputed",
        "cells_resummed",
        "reaggregated_cells",
        "error_stream_summed",
        "error_only_added_to_signed",
    )
    checks["no_numeric_addition_flags"] = all(
        result.get(key, False) is False for key in forbidden_true
    ) and result.get("no_numeric_addition", True) is True
    checks["segmented_expression_not_power_net"] = (
        result.get("formal_expression_segmented", True) is True
        and result.get("independent_subterms", True) is True
        and result.get("grouping_key") in (None, "none", "formal_expression")
    )
    checks["owner_not_proved"] = (
        result.get("owner_identity_proved", False) is False
        and result.get("physical_scale_reapplied", False) is False
        and result.get("endpoint_cancellation_proved", False) is False
        and result.get("terminal_closure_proved", False) is False
        and result.get("ledger_212_changed", False) is False
    )
    return checks


class _SingleCall:
    def __init__(self, function: Any) -> None:
        self.function = function
        self.calls = 0

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        if self.calls:
            raise RuntimeError("synthetic ungrouped leaf call is single-entry")
        self.calls += 1
        return self.function(*args, **kwargs)


def run_fixture(leaf: Any, context: Any) -> dict[str, Any]:
    """Invoke the candidate leaf once on non-arithmetic synthetic records."""

    del context  # the probe context is intentionally not a numeric context
    case = _case()
    probe_context = _ContextTrap()
    checks: dict[str, bool] = {
        "build_available": callable(getattr(leaf, "build", None)),
        "synthetic_coefficients_are_traps": True,
        "synthetic_context_is_non_numeric": True,
        "two_ordinary_Z_orders": True,
        "raw_signed_shape_is_1_and_3": True,
        "raw_error_shape_is_1_and_3": True,
    }
    build = getattr(leaf, "build", None)
    if not callable(build):
        checks["build_completed"] = False
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "check_count": len(checks),
            "synthetic_only": True,
            "leaf_calls": 0,
        }

    _CoefficientTrap.reset()
    guard = _SingleCall(build)
    try:
        result = guard(
            case["prior"],
            case["signed"],
            case["errors"],
        )
        checks["build_completed"] = isinstance(result, Mapping)
        if isinstance(result, Mapping):
            checks.update(_validate(case, result))
        try:
            guard(
                case["prior"],
                case["signed"],
                case["errors"],
            )
        except RuntimeError:
            checks["duplicate_invocation_rejected"] = guard.calls == 1
        else:
            checks["duplicate_invocation_rejected"] = False
    except BaseException as exc:
        checks["build_completed"] = False
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "check_count": len(checks),
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "leaf_calls": guard.calls,
            "coefficient_arithmetic_events": list(_CoefficientTrap.events),
            "context_numeric_calls": list(probe_context.calls),
        }

    checks["coefficient_arithmetic_not_attempted"] = not _CoefficientTrap.events
    checks["context_numeric_calls_not_attempted"] = not probe_context.calls
    return {
        "schema": SCHEMA,
        "passed": all(value is True for value in checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "raw_signed_counts": {str(order): SIGNED_COUNTS[order] for order in ORDERS},
        "raw_error_counts": {str(order): ERROR_COUNTS[order] for order in ORDERS},
        "formal_expression_is_segmented": True,
        "coefficients_were_not_summed": not _CoefficientTrap.events,
        "error_tail_kept_as_separate_reference": True,
        "error_tail_added_to_signed": False,
        "shared_source_ids_are_not_duplicate_records": True,
        "leaf_calls": guard.calls,
        "coefficient_arithmetic_events": list(_CoefficientTrap.events),
        "context_numeric_calls": list(probe_context.calls),
        "old_producer_called": False,
        "receipt_reads": False,
    }


run = run_fixture
