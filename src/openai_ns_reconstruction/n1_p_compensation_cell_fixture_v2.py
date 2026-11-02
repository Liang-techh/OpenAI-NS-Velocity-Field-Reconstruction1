"""Independent synthetic checks for the P-only physical m15 cell sum.

The fixture deliberately owns a tiny packet rather than reading a producer or
receipt.  It exercises the complete grouping key used by the companion leaf:
branch, ordinary Z order, target cell, derivative pair, formal powers,
uncertainty, and ordered error classes.  The packet contains a same-key pair
whose directed interval sum is exact, two different derivative pairs that may
not merge, an error-only record, an ordinary Z1 record, and equal records in
two different cells/branches.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any, Callable, Iterator


SCHEMA = "n1_p_compensation_cell_sum_fixture_v2"
BRANCHES = ("G0_P_W0", "G0_P_DeltaW0", "D_P_W0", "D_P_DeltaW0")
ORDERS = (0, 1)
CELLS = tuple(range(8))
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")

__all__ = ["SCHEMA", "BRANCHES", "ORDERS", "CELLS", "run_fixture", "run"]


def _q(context: Any, value: Fraction | int) -> Any:
    value = Fraction(value)
    return context.mpf(value.numerator) / value.denominator


def _box(context: Any, lower: Fraction | int, upper: Fraction | int | None = None) -> Any:
    if upper is None:
        upper = lower
    lo = _q(context, lower)
    hi = _q(context, upper)
    return context.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    context: Any,
    source_id: str,
    coefficient: Fraction | int,
    *,
    powers: tuple[int, int, int, int],
    cell: int,
    lead_order: int,
    word_order: int,
    uncertainty: bool = False,
    error_only: bool = False,
    error_classes: tuple[str, ...] = ("fixture:signed",),
) -> dict[str, Any]:
    classes = list(error_classes)
    origin: dict[str, Any] = {
        "source_ids": [source_id],
        "error_classes": classes,
        "uncertainty": uncertainty,
        "target_cell": cell,
        "lead_order": lead_order,
        "word_order": word_order,
        "symbolic_cancellation_claimed": False,
    }
    if error_only:
        origin["error_only"] = True
    result: dict[str, Any] = {
        "powers": list(powers),
        "coefficient": _box(context, coefficient),
        "origin": origin,
        "source_ids": [source_id],
        "error_classes": classes,
        "uncertainty": uncertainty,
        "target_cell": cell,
        "lead_order": lead_order,
        "word_order": word_order,
        "symbolic_cancellation_claimed": False,
    }
    if error_only:
        result["error_only"] = True
        result["symmetric_error"] = True
    return result


def _packet(context: Any) -> dict[str, dict[str, dict[str, list[dict[str, Any]]]]]:
    packet: dict[str, dict[str, dict[str, list[dict[str, Any]]]]] = {}
    for branch in BRANCHES:
        packet[branch] = {str(order): {} for order in ORDERS}
        # The same structural key is present in every branch and in cell 0;
        # branch/cell are part of the key, so no cross-branch/cell merge is
        # admissible.
        for order in ORDERS:
            for cell in CELLS:
                rows: list[dict[str, Any]] = [
                    _term(
                        context,
                        f"fixture:{branch}:cell{cell}:Z{order}:same:a",
                        Fraction(3, 16),
                        powers=(1, 0, 0, -1),
                        cell=cell,
                        lead_order=0,
                        word_order=order,
                    ),
                    _term(
                        context,
                        f"fixture:{branch}:cell{cell}:Z{order}:same:b",
                        Fraction(1, 8),
                        powers=(1, 0, 0, -1),
                        cell=cell,
                        lead_order=0,
                        word_order=order,
                    ),
                    _term(
                        context,
                        f"fixture:{branch}:cell{cell}:Z{order}:signed-uncertain",
                        Fraction(1, 8),
                        powers=(1, 0, 0, -1),
                        cell=cell,
                        lead_order=0,
                        word_order=order,
                        uncertainty=True,
                        error_classes=("fixture:signed-uncertain",),
                    ),
                    _term(
                        context,
                        f"fixture:{branch}:cell{cell}:Z{order}:error",
                        0,
                        powers=(1, 0, 0, -1),
                        cell=cell,
                        lead_order=0,
                        word_order=order,
                        uncertainty=True,
                        error_only=True,
                        error_classes=("fixture:error-only",),
                    ),
                ]
                if order == 1:
                    # Same powers and classes, but different derivative pairs.
                    rows.extend(
                        [
                            _term(
                                context,
                                f"fixture:{branch}:cell{cell}:Z1:deriv-lead",
                                Fraction(1, 5),
                                powers=(1, 0, 0, -1),
                                cell=cell,
                                lead_order=1,
                                word_order=0,
                            ),
                            _term(
                                context,
                                f"fixture:{branch}:cell{cell}:Z1:deriv-word",
                                Fraction(-1, 7),
                                powers=(1, 0, 0, -1),
                                cell=cell,
                                lead_order=0,
                                word_order=1,
                            ),
                        ]
                    )
                    rows.append(
                        _term(
                            context,
                            f"fixture:{branch}:cell{cell}:Z1:ordinary",
                            Fraction(2, 9),
                            powers=(0, 1, 0, -2),
                            cell=cell,
                            lead_order=0,
                            word_order=1,
                        )
                    )
                packet[branch][str(order)][str(cell)] = rows
    return packet


def _callable(leaf: Any) -> Callable[..., Any]:
    if callable(leaf):
        return leaf
    for name in ("build", "combine", "sum_cells", "aggregate"):
        candidate = getattr(leaf, name, None)
        if callable(candidate):
            return candidate
    raise TypeError("leaf must expose one P-cell sum callable")


def _walk(value: Any) -> Iterator[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        yield value
        for child in value.values():
            yield from _walk(child)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for child in value:
            yield from _walk(child)


def _ids(value: Mapping[str, Any]) -> set[str]:
    found: set[str] = set()
    for holder in (value, value.get("origin")):
        if not isinstance(holder, Mapping):
            continue
        raw = holder.get("source_ids", holder.get("source_id", ()))
        if isinstance(raw, str):
            found.add(raw)
        elif isinstance(raw, Sequence):
            found.update(str(item) for item in raw)
    members = value.get("members")
    if isinstance(members, Sequence) and not isinstance(members, (str, bytes, bytearray)):
        for member in members:
            if isinstance(member, Mapping):
                found.update(_ids(member))
    return found


def _coefficient(value: Mapping[str, Any]) -> Any:
    return value.get("coefficient", value.get("coefficient_interval"))


def _contains(context: Any, interval: Any, expected: Fraction | int) -> bool:
    try:
        point = _q(context, expected)
        return bool(interval.a <= point.a and point.b <= interval.b)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _candidate_values(result: Mapping[str, Any], ids: set[str]) -> list[Mapping[str, Any]]:
    values: list[Mapping[str, Any]] = []
    for item in _walk(result):
        if ids.issubset(_ids(item)):
            values.append(item)
    return values


def run_fixture(leaf: Any, context: Any) -> dict[str, Any]:
    """Call the candidate once on the synthetic packet and validate grouping."""

    packet = _packet(context)
    checks: dict[str, bool] = {
        "callable_available": False,
        "call_completed": False,
        "result_mapping": False,
        "same_key_sum": False,
        "different_derivative_pairs_not_merged": False,
        "signed_error_separation": False,
        "ordinary_Z1_retained": False,
        "cell_branch_separation": False,
        "formal_powers_retained": False,
        "source_metadata_retained": False,
        "no_cancellation_claim": False,
    }
    try:
        call = _callable(leaf)
        checks["callable_available"] = True
        # The minimal two positional arguments are intentional: the actual
        # runner supplies provenance/physical metadata through keyword-only
        # arguments, while a fixture must remain independent of saved inputs.
        result = call(context, packet)
        checks["call_completed"] = True
    except BaseException as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "check_count": len(checks),
            "synthetic_only": True,
            "error": f"{type(exc).__name__}: {exc}",
        }
    checks["result_mapping"] = isinstance(result, Mapping)
    if not checks["result_mapping"]:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "check_count": len(checks),
            "synthetic_only": True,
        }
    result_map = result
    all_values = list(_walk(result_map))
    all_terms = [item for item in all_values if _coefficient(item) is not None]
    checks["formal_powers_retained"] = all(
        isinstance(item.get("powers"), Sequence)
        and len(item["powers"]) == 4
        and all(type(power) is int for power in item["powers"])
        for item in all_terms
    )
    checks["source_metadata_retained"] = any(
        {"source_ids", "error_classes", "uncertainty"}.issubset(item.get("origin", {}))
        for item in all_terms
        if isinstance(item.get("origin"), Mapping)
    )

    same_ids = {
        f"fixture:G0_P_W0:cell0:Z0:same:a",
        f"fixture:G0_P_W0:cell0:Z0:same:b",
    }
    same_candidates = _candidate_values(result_map, same_ids)
    checks["same_key_sum"] = any(
        _coefficient(item) is not None and _contains(context, _coefficient(item), Fraction(5, 16))
        for item in same_candidates
    )

    derivative_ids = {
        "fixture:G0_P_W0:cell0:Z1:deriv-lead",
        "fixture:G0_P_W0:cell0:Z1:deriv-word",
    }
    checks["different_derivative_pairs_not_merged"] = not any(
        derivative_ids.issubset(_ids(item)) for item in all_values
    )
    # A signed aggregate must not contain the error-only source alongside a
    # signed source.  The error source itself may be retained in a separate
    # stream or member reference.
    error_id = "fixture:G0_P_W0:cell0:Z0:error"
    signed_id = "fixture:G0_P_W0:cell0:Z0:same:a"
    uncertain_signed_id = "fixture:G0_P_W0:cell0:Z0:signed-uncertain"
    checks["signed_error_separation"] = not any(
        error_id in _ids(item) and signed_id in _ids(item)
        for item in all_values
    ) and any(uncertain_signed_id in _ids(item) for item in all_values)
    checks["ordinary_Z1_retained"] = any(
        any(source.endswith(":Z1:ordinary") for source in _ids(item))
        for item in all_values
    )
    cell0_id = "fixture:G0_P_W0:cell0:Z0:same:a"
    cell1_id = "fixture:G0_P_W0:cell1:Z0:same:a"
    branch2_id = "fixture:D_P_W0:cell0:Z0:same:a"
    checks["cell_branch_separation"] = not any(
        (cell0_id in _ids(item) and cell1_id in _ids(item))
        or (cell0_id in _ids(item) and branch2_id in _ids(item))
        for item in all_values
    )
    checks["no_cancellation_claim"] = all(
        item.get("symbolic_cancellation_claimed", False) is False
        and item.get("cancellation_performed", False) is False
        for item in all_values
    )
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "exact_rational_inputs": True,
        "same_key_sum_expected": "3/16 + 1/8 = 5/16",
        "error_only_is_not_signed": True,
        "ordinary_Z1_checked": True,
        "cell_and_branch_keys_checked": True,
        "cross_cell_cancellation": False,
        "cross_branch_cancellation": False,
    }


def run(leaf: Any, context: Any) -> dict[str, Any]:
    return run_fixture(leaf, context)

