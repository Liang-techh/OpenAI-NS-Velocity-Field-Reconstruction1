"""Independent exact-rational fixture for first-collar sector grouping.

The fixture supplies already-weighted synthetic sectors to the new leaf.  It
does not read a receipt or import a producer.  Singleton rational values and
finite rational interval boxes exercise first-seen grouping, ordinary Z1,
provenance retention, sign retention, and the separation of unequal error
keys.  An overlapping interval is deliberately used only as an interval
example; it is never treated as evidence that two sources are equal.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "first_collar_m15_group_fixture_v1"
SOURCE_ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")

__all__ = ["SCHEMA", "SOURCE_ATOMS", "run_fixture", "run"]


def _q(context: Any, value: Fraction | int) -> Any:
    rational = Fraction(value)
    return context.mpf(rational.numerator) / rational.denominator


def _interval(
    context: Any,
    lower: Fraction | int,
    upper: Fraction | int | None = None,
) -> Any:
    if upper is None:
        upper = lower
    lo = _q(context, lower)
    hi = _q(context, upper)
    lo_raw = getattr(lo, "_mpi_", None)
    hi_raw = getattr(hi, "_mpi_", None)
    if lo_raw is not None and hi_raw is not None and hasattr(context, "make_mpf"):
        return context.make_mpf((lo_raw[0], hi_raw[1]))
    return context.mpf([lo, hi])


def _contains_float(value: Any) -> bool:
    if isinstance(value, float):
        return True
    if isinstance(value, Mapping):
        return any(_contains_float(key) or _contains_float(item) for key, item in value.items())
    if isinstance(value, (list, tuple, set, frozenset)):
        return any(_contains_float(item) for item in value)
    return False


def _inside(context: Any, value: Any, expected: Fraction | int) -> bool:
    try:
        point = _q(context, expected)
        return bool(value.a <= point.a and point.b <= value.b)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _positive(value: Any) -> bool:
    try:
        return bool(value.a > 0)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _negative(value: Any) -> bool:
    try:
        return bool(value.b < 0)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _overlap(left: Any, right: Any) -> bool:
    try:
        return bool(left.a <= right.b and right.a <= left.b)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _term(
    context: Any,
    *,
    source_id: str,
    powers: tuple[int, int, int, int],
    lower: Fraction | int,
    upper: Fraction | int | None = None,
    uncertainty: bool,
    error_class: str,
) -> dict[str, Any]:
    return {
        "powers": list(powers),
        "coefficient": _interval(context, lower, upper),
        "origin": {
            "source_ids": [source_id],
            "uncertainty": uncertainty,
            "error_classes": [error_class],
        },
    }


def _packet(context: Any) -> dict[str, dict[int, list[dict[str, Any]]]]:
    # The first two Z0 records have overlapping boxes and the same structural
    # key.  They must merge with the third same-key record from swirl.
    z0_axial = [
        _term(
            context,
            source_id="fixture:z0:axial:singleton",
            powers=(1, 0, 0, 0),
            lower=Fraction(1, 3),
            uncertainty=False,
            error_class="shared-leading",
        ),
        _term(
            context,
            source_id="fixture:z0:axial:overlap",
            powers=(1, 0, 0, 0),
            lower=Fraction(1, 4),
            upper=Fraction(1, 2),
            uncertainty=False,
            error_class="shared-leading",
        ),
    ]
    z0_swirl = [
        _term(
            context,
            source_id="fixture:z0:swirl:shared-key",
            powers=(1, 0, 0, 0),
            lower=Fraction(1, 2),
            uncertainty=False,
            error_class="shared-leading",
        )
    ]
    z0_pressure = [
        _term(
            context,
            source_id="fixture:z0:pressure:regular-error",
            powers=(1, 0, 0, 0),
            lower=Fraction(-1, 7),
            uncertainty=False,
            error_class="pressure-error",
        ),
        _term(
            context,
            source_id="fixture:z0:pressure:uncertain-error",
            powers=(1, 0, 0, 0),
            lower=Fraction(1, 9),
            uncertainty=True,
            error_class="pressure-error-uncertain",
        ),
    ]

    z1_axial = [
        _term(
            context,
            source_id="fixture:z1:axial:shared-key",
            powers=(2, 1, 0, -1),
            lower=Fraction(-1, 3),
            uncertainty=False,
            error_class="ordinary-z1",
        )
    ]
    z1_swirl = [
        _term(
            context,
            source_id="fixture:z1:swirl:shared-key",
            powers=(2, 1, 0, -1),
            lower=Fraction(-1, 6),
            uncertainty=False,
            error_class="ordinary-z1",
        )
    ]
    z1_pressure = [
        _term(
            context,
            source_id="fixture:z1:pressure:shared-key",
            powers=(2, 1, 0, -1),
            lower=Fraction(1, 10),
            uncertainty=False,
            error_class="ordinary-z1",
        ),
        _term(
            context,
            source_id="fixture:z1:pressure:uncertain-key",
            powers=(2, 1, 0, -1),
            lower=Fraction(1, 100),
            uncertainty=True,
            error_class="ordinary-z1-uncertain",
        ),
    ]
    return {
        "axial": {0: z0_axial, 1: z1_axial},
        "swirl": {0: z0_swirl, 1: z1_swirl},
        "pressure": {0: z0_pressure, 1: z1_pressure},
    }


def _all_groups(result: Mapping[str, Any], order: int) -> list[Mapping[str, Any]]:
    groups = result.get("groups")
    if not isinstance(groups, Mapping):
        return []
    row = groups.get(order, groups.get(str(order), ()))
    return list(row) if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)) else []


def _find_group(
    groups: Sequence[Mapping[str, Any]],
    *,
    uncertainty: bool,
    error_class: str,
) -> Mapping[str, Any] | None:
    for group in groups:
        origin = group.get("origin")
        if not isinstance(origin, Mapping):
            continue
        if origin.get("uncertainty") is not uncertainty:
            continue
        if origin.get("error_classes") == [error_class]:
            return group
    return None


def _member_keys(result: Mapping[str, Any]) -> set[tuple[Any, ...]]:
    keys: set[tuple[Any, ...]] = set()
    for order in (0, 1):
        for group in _all_groups(result, order):
            members = group.get("members", ())
            if not isinstance(members, Sequence):
                continue
            for member in members:
                if isinstance(member, Mapping):
                    keys.add((member.get("channel"), member.get("ordinary_Z_order"), member.get("index")))
    return keys


def run_fixture(leaf: Any, context: Any) -> dict[str, Any]:
    """Run one leaf call against exact synthetic sectors and inspect it."""

    combine = getattr(leaf, "combine", None)
    if not callable(combine):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"leaf_combine_available": False},
            "check_count": 1,
            "synthetic_only": True,
            "implementation_imported": False,
        }

    packet = _packet(context)
    checks: dict[str, bool] = {
        "exact_packet_has_no_float": not _contains_float(packet),
        "source_atoms_are_canonical": SOURCE_ATOMS == (
            "hb",
            "F0_anchor",
            "Pstar_squared",
            "Lambda",
        ),
    }
    try:
        result = combine(
            packet["axial"],
            packet["swirl"],
            packet["pressure"],
            source_atoms=SOURCE_ATOMS,
        )
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {**checks, "combine_completed": False},
            "check_count": len(checks) + 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["combine_completed"] = isinstance(result, Mapping)
    checks["basis_retained"] = result.get("source_power_basis") == list(SOURCE_ATOMS)
    checks["ordinary_Z0_Z1_retained"] = result.get("ordinary_Z_orders") == [0, 1]
    checks["all_channels_retained"] = result.get("channels") == ["axial", "swirl", "pressure"]
    checks["physical_factors_not_applied"] = result.get("physical_factors_applied") is False
    checks["cancellation_not_performed"] = result.get("cancellation_performed") is False
    checks["interval_addition_only"] = result.get("outward_interval_addition_only") is True

    z0 = _all_groups(result, 0)
    z1 = _all_groups(result, 1)
    checks["multiple_same_key_records_merge"] = len(z0) == 3 and len(z1) == 2

    shared0 = _find_group(z0, uncertainty=False, error_class="shared-leading")
    shared1 = _find_group(z1, uncertainty=False, error_class="ordinary-z1")
    regular0 = _find_group(z0, uncertainty=False, error_class="pressure-error")
    uncertain0 = _find_group(z0, uncertainty=True, error_class="pressure-error-uncertain")
    uncertain1 = _find_group(z1, uncertainty=True, error_class="ordinary-z1-uncertain")
    checks["Z0_shared_group_has_three_members"] = (
        shared0 is not None and shared0.get("origin", {}).get("member_count") == 3
        and len(shared0.get("members", ())) == 3
    )
    checks["Z1_shared_group_has_three_members"] = (
        shared1 is not None and shared1.get("origin", {}).get("member_count") == 3
        and len(shared1.get("members", ())) == 3
    )
    checks["unequal_error_keys_are_separate"] = (
        regular0 is not None and uncertain0 is not None and regular0 is not uncertain0
        and uncertain1 is not None and shared1 is not uncertain1
    )
    checks["Z0_shared_sign_and_value"] = (
        shared0 is not None and _positive(shared0["coefficient"])
        and _inside(context, shared0["coefficient"], Fraction(5, 4))
    )
    checks["Z0_negative_sign_retained"] = regular0 is not None and _negative(regular0["coefficient"])
    checks["Z1_shared_nonzero_negative_sign"] = (
        shared1 is not None and _negative(shared1["coefficient"])
        and _inside(context, shared1["coefficient"], Fraction(-2, 5))
    )
    checks["Z1_uncertain_positive_sign"] = uncertain1 is not None and _positive(uncertain1["coefficient"])

    expected_member_keys = {
        (channel, order, index)
        for channel, rows in packet.items()
        for order, terms in rows.items()
        for index, _ in enumerate(terms)
    }
    checks["every_input_member_retained_once"] = (
        _member_keys(result) == expected_member_keys
        and sum(
            len(group.get("members", ()))
            for order in (0, 1)
            for group in _all_groups(result, order)
        ) == len(expected_member_keys)
    )

    provenance_ok = True
    for order in (0, 1):
        for group in _all_groups(result, order):
            origin = group.get("origin", {})
            if not isinstance(origin, Mapping):
                provenance_ok = False
                continue
            for member in group.get("members", ()):
                if not isinstance(member, Mapping):
                    provenance_ok = False
                    continue
                channel = member.get("channel")
                index = member.get("index")
                original = packet.get(channel, {}).get(order, [])[index] if channel in packet else None
                provenance = member.get("provenance")
                if not isinstance(original, Mapping) or not isinstance(provenance, Mapping):
                    provenance_ok = False
                    continue
                provenance_ok = provenance_ok and (
                    provenance.get("source_ids") == original["origin"]["source_ids"]
                    and provenance.get("uncertainty") is original["origin"]["uncertainty"]
                    and provenance.get("error_classes") == original["origin"]["error_classes"]
                    and member.get("record", {}).get("origin") == original["origin"]
                )
    checks["source_ids_errors_and_records_retained"] = provenance_ok
    checks["ordinary_Z1_members_are_labeled"] = all(
        member.get("ordinary_Z_order") == 1
        for group in z1
        for member in group.get("members", ())
    )

    overlap_seen = _overlap(
        packet["axial"][0][0]["coefficient"],
        packet["axial"][0][1]["coefficient"],
    )
    checks["overlap_is_exercised"] = overlap_seen
    checks["overlap_does_not_prove_source_equality"] = (
        "source_equality_proved" not in result
        and "source_equality_from_overlap" not in result
    )

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "exact_rational_inputs": True,
        "exact_rational_singletons_and_boxes": True,
        "mpmath_directed_intervals": True,
        "ordinary_Z1_checked": True,
        "member_records_retained": True,
        "provenance_retained": True,
        "no_float_inputs": True,
        "source_equality_proved": False,
        "source_equality_from_overlap": False,
        "implementation_imported": False,
        "old_code_called": False,
    }


def run(leaf: Any, context: Any) -> dict[str, Any]:
    return run_fixture(leaf, context)

