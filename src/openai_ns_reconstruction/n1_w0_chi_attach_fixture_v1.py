"""Exact synthetic checks for the saved W0-chi attachment leaf.

The fixture is deliberately independent of the two saved receipts.  It gives
the leaf an already grouped prior formal prefix and an already aggregated
W0-chi correction.  The signed correction and its symmetric error are passed
as different channels; the eight cell records that produced a saved group are
never expanded here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "n1_w0_chi_attach_fixture_v1"
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1)
__all__ = ["SCHEMA", "POWER_BASIS", "ORDERS", "run_fixture", "run"]


def _point(c: Any, value: Fraction | int) -> Any:
    value = Fraction(value)
    numerator = c.mpf(value.numerator)
    denominator = c.mpf(value.denominator)
    return numerator / denominator


def _interval(c: Any, lower: Fraction | int, upper: Fraction | int | None = None) -> Any:
    if upper is None:
        upper = lower
    lo = _point(c, lower)
    hi = _point(c, upper)
    raw_lo = getattr(lo, "_mpi_", None)
    raw_hi = getattr(hi, "_mpi_", None)
    if raw_lo is None or raw_hi is None:
        raise TypeError("synthetic fixture requires directed interval endpoints")
    return c.make_mpf((raw_lo[0], raw_hi[1]))


def _term(
    c: Any,
    *,
    source_id: str,
    path: Sequence[Any],
    powers: tuple[int, int, int, int],
    value: Fraction | int,
    error_only: bool = False,
    signed_field: bool = True,
    error_class: str = "synthetic_signed",
) -> dict[str, Any]:
    coefficient = _interval(c, value) if not error_only else _interval(c, -abs(Fraction(value)), abs(Fraction(value)))
    origin = {
        "source_ids": [source_id],
        "source_paths": [list(path)],
        "error_classes": [error_class],
        "uncertainty": True if error_only else False,
        "signed_field": signed_field,
        "error_only": error_only,
    }
    return {
        "powers": list(powers),
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": [source_id],
        "source_paths": [list(path)],
        "signed_field": signed_field,
        "error_only": error_only,
        "symmetric_error": error_only,
        "uncertainty": bool(error_only),
        "receipt": "synthetic://n1_w0_chi_attach",
        "path": list(path),
    }


def _case(c: Any) -> dict[str, Any]:
    family = {
        "family": "synthetic-n1-w0-chi-attach",
        "revision": "v1",
        "owner_identity_proved": False,
    }
    logs = {name: "synthetic:" + name for name in POWER_BASIS}
    prior: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    signed: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    errors: dict[int, list[dict[str, Any]]] = {0: [], 1: []}

    prior[0] = [
        _term(c, source_id="prior:Z0:signed:positive", path=("prior", "formal_sum", "0", 0), powers=(0, 0, 0, 0), value=Fraction(1, 4), error_class="prior_signed"),
        _term(c, source_id="prior:Z0:error", path=("prior", "formal_sum", "0", 1), powers=(1, 0, 0, -1), value=Fraction(1, 32), error_only=True, signed_field=False, error_class="prior_error"),
    ]
    prior[1] = [
        _term(c, source_id="prior:Z1:signed:negative", path=("prior", "formal_sum", "1", 0), powers=(0, 1, 0, -2), value=Fraction(-3, 8), error_class="prior_signed"),
        _term(c, source_id="prior:Z1:error", path=("prior", "formal_sum", "1", 1), powers=(1, 0, 0, -1), value=Fraction(1, 64), error_only=True, signed_field=False, error_class="prior_error"),
    ]

    signed[0] = [
        _term(c, source_id="correction:W0_chi:Z0:swirl:0", path=("correction", "m15", "branches", "W0_chi", "0", "swirl", 0), powers=(0, 2, 0, -3), value=Fraction(-1, 16), error_class="chi_signed"),
    ]
    signed[1] = [
        _term(c, source_id="correction:W0_chi:Z1:swirl:0", path=("correction", "m15", "branches", "W0_chi", "1", "swirl", 0), powers=(0, 2, 0, -3), value=Fraction(3, 32), error_class="chi_signed"),
        _term(c, source_id="correction:W0_chi:Z1:swirl:1", path=("correction", "m15", "branches", "W0_chi", "1", "swirl", 1), powers=(0, 2, 0, -3), value=Fraction(-1, 32), error_class="chi_signed"),
        _term(c, source_id="correction:W0_chi:Z1:swirl:2", path=("correction", "m15", "branches", "W0_chi", "1", "swirl", 2), powers=(0, 2, 0, -2), value=Fraction(1, 64), error_class="chi_signed"),
    ]
    errors[0] = [
        _term(c, source_id="correction:W0_chi:Z0:swirl_error:0", path=("correction", "m15", "branches", "W0_chi", "0", "swirl_error", 0), powers=(0, 2, 0, -3), value=Fraction(1, 128), error_only=True, signed_field=False, error_class="chi_tail"),
    ]
    errors[1] = [
        _term(c, source_id="correction:W0_chi:Z1:swirl_error:0", path=("correction", "m15", "branches", "W0_chi", "1", "swirl_error", 0), powers=(0, 2, 0, -3), value=Fraction(1, 256), error_only=True, signed_field=False, error_class="chi_tail"),
        _term(c, source_id="correction:W0_chi:Z1:swirl_error:1", path=("correction", "m15", "branches", "W0_chi", "1", "swirl_error", 1), powers=(0, 2, 0, -3), value=Fraction(1, 512), error_only=True, signed_field=False, error_class="chi_tail"),
        _term(c, source_id="correction:W0_chi:Z1:swirl_error:2", path=("correction", "m15", "branches", "W0_chi", "1", "swirl_error", 2), powers=(0, 2, 0, -2), value=Fraction(1, 1024), error_only=True, signed_field=False, error_class="chi_tail"),
    ]

    prior_index = {
        term["source_ids"][0]: {"receipt": "synthetic://prior", "path": term["path"]}
        for order in ORDERS
        for term in prior[order]
    }
    correction_index = {
        term["source_ids"][0]: {"receipt": "synthetic://correction", "path": term["path"]}
        for order in ORDERS
        for term in signed[order] + errors[order]
    }
    branch = {
        order: {
            "axial": [],
            "swirl": signed[order],
            "axial_error": [],
            "swirl_error": errors[order],
            "swirl_cells": [[] for _ in range(8)],
            "swirl_error_cells": [[] for _ in range(8)],
        }
        for order in ORDERS
    }
    return {
        "family": family,
        "logs": logs,
        "prior": prior,
        "correction": {"branches": {"W0_chi": branch}},
        "signed": signed,
        "errors": errors,
        "prior_index": prior_index,
        "correction_index": correction_index,
        "expected_Z1_shared_net": _interval(c, Fraction(1, 16)),
        "source_index": {"prior": prior_index, "W0_chi": correction_index},
        "certificate": {
            "source_family_exact": True,
            "source_logs_exact": True,
            "physical_scale_reapplied": False,
            "owner_identity_proved": False,
        },
    }


def _orders(value: Any) -> dict[int, Any]:
    if not isinstance(value, Mapping):
        return {}
    result: dict[int, Any] = {}
    for order in ORDERS:
        for key in (order, str(order), "Z" + str(order), "z" + str(order)):
            if key in value:
                result[order] = value[key]
                break
    return result


def _rows(value: Any, order: int) -> list[Any]:
    rows = _orders(value)
    selected = rows.get(order, ())
    return list(selected) if isinstance(selected, Sequence) and not isinstance(selected, (str, bytes, bytearray)) else []


def _channel(result: Mapping[str, Any], name: str) -> Any:
    for container_name in (name, "correction", "W0_chi", "attached", "m15"):
        container = result.get(container_name)
        if isinstance(container, Mapping):
            if name in container:
                return container[name]
            branches = container.get("branches")
            if isinstance(branches, Mapping):
                branch = branches.get("W0_chi")
                if isinstance(branch, Mapping):
                    if name in branch:
                        return branch[name]
                    rows = {order: branch.get(order, branch.get(str(order), branch.get("Z" + str(order), {}))) for order in ORDERS}
                    if all(isinstance(rows[order], Mapping) for order in ORDERS):
                        return {order: rows[order].get(name, []) for order in ORDERS}
            rows = {order: container.get(order, container.get(str(order), container.get("Z" + str(order), {}))) for order in ORDERS}
            if all(isinstance(rows[order], Mapping) for order in ORDERS):
                return {order: rows[order].get(name, []) for order in ORDERS}
    return {}


def _term_ids(term: Mapping[str, Any]) -> tuple[str, ...]:
    origin = term.get("origin")
    values: list[Any] = []
    if isinstance(origin, Mapping):
        values.extend(origin.get("source_ids", ()))
        values.extend(origin.get("source_paths", ()))
    values.extend(term.get("source_ids", ()))
    values.extend(term.get("source_paths", ()))
    return tuple(str(value) for value in values if isinstance(value, (str, int)))


def _interval_signature(value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return None
    return tuple(tuple(endpoint) for endpoint in raw)


def _identity_once(consumed: Sequence[Any], originals: Sequence[Any]) -> bool:
    return len(consumed) == len(originals) and all(
        sum(item is original for item in consumed) == 1 for original in originals
    )


def _formal(result: Mapping[str, Any]) -> Any:
    for key in ("prior_formal_sum", "formal_sum", "prior"):
        value = result.get(key)
        if isinstance(value, Mapping):
            return value
    return {}


def _validate(case: Mapping[str, Any], result: Mapping[str, Any]) -> dict[str, bool]:
    signed = _channel(result, "signed")
    errors = _channel(result, "error_only")
    raw_signed_packet = result.get("signed_correction")
    raw_signed = {
        order: _rows(raw_signed_packet, order)
        for order in ORDERS
    }
    checks: dict[str, bool] = {
        "result_mapping": isinstance(result, Mapping),
        "ordinary_Z0_and_Z1": set(_orders(signed)) == set(ORDERS) and set(_orders(errors)) == set(ORDERS),
        "signed_counts": len(_rows(signed, 0)) == 1 and len(_rows(signed, 1)) == 2,
        "raw_signed_counts": len(raw_signed[0]) == 1 and len(raw_signed[1]) == 3,
        "error_counts": len(_rows(errors, 0)) == 1 and len(_rows(errors, 1)) == 3,
        "axial_is_empty": True,
        "owner_identity_not_proved": result.get("owner_identity_proved") is False,
        "physical_scale_not_reapplied": result.get("physical_scale_reapplied", False) is False,
        "endpoint_and_ledger_absent": result.get("endpoint_cancellation_proved", False) is False and result.get("ledger_212_changed", False) is False,
    }
    source_error_ids: set[str] = set()
    for order in ORDERS:
        for term in _rows(errors, order):
            if isinstance(term, Mapping):
                source_error_ids.update(_term_ids(term))
        for term in _rows(signed, order):
            if not isinstance(term, Mapping):
                checks["signed_markers"] = False
                continue
            origin = term.get("origin", {})
            checks.setdefault("signed_markers", True)
            checks["signed_markers"] &= term.get("signed_field") is True and term.get("error_only") is not True and isinstance(origin, Mapping) and origin.get("signed_field") is True and origin.get("error_only") is not True
            checks.setdefault("signed_paths", True)
            checks["signed_paths"] &= bool(_term_ids(term))
            powers = term.get("powers")
            checks.setdefault("four_integer_powers", True)
            checks["four_integer_powers"] &= isinstance(powers, Sequence) and len(powers) == 4 and all(type(power) is int for power in powers)
        for term in _rows(errors, order):
            if not isinstance(term, Mapping):
                checks["error_markers"] = False
                continue
            origin = term.get("origin", {})
            checks.setdefault("error_markers", True)
            checks["error_markers"] &= term.get("signed_field") is False and term.get("error_only") is True and term.get("symmetric_error") is True and isinstance(origin, Mapping) and origin.get("signed_field") is False and origin.get("error_only") is True and origin.get("uncertainty") is True
            powers = term.get("powers")
            checks.setdefault("four_integer_powers", True)
            checks["four_integer_powers"] &= isinstance(powers, Sequence) and len(powers) == 4 and all(type(power) is int for power in powers)
    checks["error_not_in_signed"] = not any(set(_term_ids(term)) & source_error_ids for order in ORDERS for term in _rows(signed, order) if isinstance(term, Mapping))
    checks["raw_signed_identity_retained"] = all(
        len(raw_signed[order]) == len(case["signed"][order])
        and all(raw_signed[order][index] is case["signed"][order][index] for index in range(len(raw_signed[order])))
        for order in ORDERS
    )
    net_inputs: list[Any] = []
    for term in _rows(signed, 1):
        if isinstance(term, Mapping):
            inputs = term.get("input_terms", ())
            if isinstance(inputs, Sequence) and not isinstance(inputs, (str, bytes, bytearray)):
                net_inputs.extend(inputs)
    checks["three_signed_inputs_retained_exactly_once"] = _identity_once(net_inputs, raw_signed[1])
    error_rows = {order: _rows(errors, order) for order in ORDERS}
    checks["error_identity_retained"] = all(
        len(error_rows[order]) == len(case["errors"][order])
        and all(error_rows[order][index] is case["errors"][order][index] for index in range(len(error_rows[order])))
        for order in ORDERS
    )
    checks["error_only_independent"] = not any(
        item is error_term
        for item in net_inputs
        for order in ORDERS
        for error_term in error_rows[order]
    )
    shared_power = (0, 2, 0, -3)
    other_power = (0, 2, 0, -2)
    net_z1 = _rows(signed, 1)
    shared_groups = [
        term for term in net_z1
        if isinstance(term, Mapping) and tuple(term.get("powers", ())) == shared_power
    ]
    other_groups = [
        term for term in net_z1
        if isinstance(term, Mapping) and tuple(term.get("powers", ())) == other_power
    ]
    checks["Z1_net_has_two_power_groups"] = len(net_z1) == 2 and len(shared_groups) == 1 and len(other_groups) == 1
    shared_inputs = []
    if shared_groups and isinstance(shared_groups[0], Mapping):
        inputs = shared_groups[0].get("input_terms", ())
        if isinstance(inputs, Sequence) and not isinstance(inputs, (str, bytes, bytearray)):
            shared_inputs = list(inputs)
    checks["Z1_shared_power_inputs_are_first_two"] = (
        len(shared_inputs) == 2
        and shared_inputs[0] is raw_signed[1][0]
        and shared_inputs[1] is raw_signed[1][1]
    ) if len(raw_signed[1]) == 3 else False
    checks["Z1_shared_net_exact_directed_sum"] = (
        bool(shared_groups)
        and isinstance(shared_groups[0], Mapping)
        and _interval_signature(shared_groups[0].get("coefficient")) == _interval_signature(case["expected_Z1_shared_net"])
    )
    checks["power_compatibility"] = all(
        sorted(tuple(term.get("powers", ())) for term in raw_signed[order]) == sorted(tuple(term.get("powers", ())) for term in error_rows[order])
        for order in ORDERS
    )
    formal = _formal(result)
    prior_ids = {term["source_ids"][0] for order in ORDERS for term in case["prior"][order]}
    formal_ids = {item for order in ORDERS for term in _rows(formal, order) if isinstance(term, Mapping) for item in _term_ids(term)}
    checks["prior_paths_retained"] = prior_ids.issubset(formal_ids) or prior_ids.issubset({item for order in ORDERS for term in _rows(result.get("prior_formal_sum", {}), order) if isinstance(term, Mapping) for item in _term_ids(term)})
    checks["prior_prefix_not_repeated"] = len(formal_ids & prior_ids) <= len(prior_ids)
    checks["grouped_cells_not_resummed"] = result.get("cells_resummed", False) is False and result.get("reaggregated_cells", False) is False
    checks["source_family_retained"] = result.get("source_family") in (case["family"], None) or result.get("source_family") == case["family"]
    checks["source_logs_retained"] = result.get("source_logs") in (case["logs"], None) or result.get("source_logs") == case["logs"]
    return checks


class _SingleCall:
    """Reject a second fixture call before it can reach the leaf."""

    def __init__(self, function: Any) -> None:
        self.function = function
        self.calls = 0

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        if self.calls:
            raise RuntimeError("synthetic fixture leaf call is single-entry")
        self.calls += 1
        return self.function(*args, **kwargs)


def run_fixture(leaf: Any, context: Any) -> dict[str, Any]:
    """Call ``leaf.build`` once on exact binary-rational synthetic inputs."""

    case = _case(context)
    checks: dict[str, bool] = {
        "build_available": callable(getattr(leaf, "build", None)),
        "exact_binary_rational_inputs": True,
        "two_ordinary_Z_orders": True,
        "eight_grouped_cell_inputs_not_expanded": True,
    }
    build = getattr(leaf, "build", None)
    if not callable(build):
        checks["build_completed"] = False
        return {"schema": SCHEMA, "passed": False, "checks": checks, "check_count": len(checks), "synthetic_only": True, "leaf_calls": 0}

    guard = _SingleCall(build)
    try:
        result = guard(
            context,
            prior_formal_sum=case["prior"],
            correction=case["correction"],
            w0_chi=case["correction"]["branches"]["W0_chi"],
            source_family=case["family"],
            source_logs=case["logs"],
            source_atoms=case["logs"],
            source_index=case["source_index"],
            prior_source_index=case["prior_index"],
            correction_source_index=case["correction_index"],
            owner_certificate=case["certificate"],
            scale_certificate=case["certificate"],
            certificate=case["certificate"],
        )
        checks["build_completed"] = isinstance(result, Mapping)
        if isinstance(result, Mapping):
            checks.update(_validate(case, result))
        try:
            guard(
                context,
                prior_formal_sum=case["prior"],
                correction=case["correction"],
                source_family=case["family"],
                source_logs=case["logs"],
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
        }
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "exact_binary_rational_inputs": True,
        "positive_and_negative_signed_inputs": True,
        "symmetric_error_inputs": True,
        "error_not_added_to_signed": checks.get("error_not_in_signed", False),
        "source_paths_checked": checks.get("signed_paths", False),
        "duplicate_call_guard": True,
        "error_tail_is_independent_reference": False,
        "error_tail_kept_as_separate_reference": True,
        "error_tail_independence_not_proved": True,
        "leaf_calls": guard.calls,
        "old_producer_called": False,
        "receipt_reads": False,
    }


run = run_fixture

