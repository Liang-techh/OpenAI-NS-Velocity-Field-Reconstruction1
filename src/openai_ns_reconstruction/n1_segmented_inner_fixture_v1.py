"""Small independent fixture for the segmented-inner candidate.

The fixture is synthetic only. Core rows are supplied directly in ordinary
``Z0``/``Z1`` order rows. Each collar input row is one already-grouped prior
record with its saved coefficient and its complete prior member list nested
inside the record; prior members are never expanded into new leaf rows.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "n1_segmented_inner_fixture_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1)

__all__ = ["SCHEMA", "SOURCE_POWER_BASIS", "run_fixture", "run"]


def _point(context: Any, value: Fraction | int) -> Any:
    fraction = Fraction(value)
    return context.mpf(fraction.numerator) / fraction.denominator


def _interval(context: Any, value: Fraction | int) -> Any:
    point = _point(context, value)
    raw = getattr(point, "_mpi_", None)
    if raw is not None and hasattr(context, "make_mpf"):
        return context.make_mpf((raw[0], raw[1]))
    return context.mpf([point, point])


def _same_interval(left: Any, right: Any) -> bool:
    return getattr(left, "_mpi_", None) == getattr(right, "_mpi_", None)


def _term(
    context: Any,
    *,
    source_id: str,
    powers: tuple[int, int, int, int],
    value: Fraction | int,
    uncertainty: bool,
    error_classes: Sequence[str],
    temporal_tag: str,
    channel: str,
    order: int,
    cell: int | None,
    word: str,
    receipt: str,
    path: list[Any],
) -> dict[str, Any]:
    origin = {
        "source_ids": [source_id],
        "uncertainty": uncertainty,
        "error_classes": list(error_classes),
    }
    return {
        "powers": list(powers),
        "F0_power": powers[1],
        "Lambda_power": powers[3],
        "coefficient": _interval(context, value),
        "origin": origin,
        "temporal_tag": temporal_tag,
        "channel": channel,
        "ordinary_Z_order": order,
        "cell": cell,
        "word": word,
        "segment": "collar" if word.startswith("collar") else "core",
        "receipt": receipt,
        "path": path,
    }


def _prior_group(
    context: Any,
    *,
    order: int,
    group_index: int,
    powers: tuple[int, int, int, int],
    group_value: Fraction,
    member_values: tuple[Fraction, Fraction],
    uncertainty: bool,
    error_classes: tuple[str, ...],
) -> dict[str, Any]:
    members: list[dict[str, Any]] = []
    source_ids: list[str] = []
    for member_index, member_value in enumerate(member_values):
        source_id = f"fixture:collar:Z{order}:group{group_index}:member{member_index}"
        source_ids.append(source_id)
        member_path = [
            "synthetic://collar", "groups", str(order), group_index,
            "members", member_index, "record",
        ]
        member_record = _term(
            context,
            source_id=source_id,
            powers=powers,
            value=member_value,
            uncertainty=uncertainty,
            error_classes=error_classes,
            temporal_tag=f"collar:prior:Z{order}:group{group_index}:member{member_index}",
            channel="collar_group",
            order=order,
            cell=None,
            word="collar-prior-member",
            receipt="synthetic://collar",
            path=member_path,
        )
        members.append({
            "channel": "collar_group",
            "ordinary_Z_order": order,
            "index": member_index,
            "record": member_record,
            "origin": member_record["origin"],
            "provenance": member_record["origin"],
            "receipt": "synthetic://collar",
            "path": member_path,
        })
    group_path = ["synthetic://collar", "groups", str(order), group_index]
    group_origin = {
        "source_ids": source_ids,
        "uncertainty": uncertainty,
        "error_classes": list(error_classes),
        "member_count": len(members),
        "provenance_retained": True,
        "source_equality_proved": False,
    }
    return {
        "group_id": f"fixture:collar:Z{order}:group{group_index}",
        "powers": list(powers),
        "coefficient": _interval(context, group_value),
        "origin": group_origin,
        "prior_group_members": members,
        "prior_group_key": {
            "powers": list(powers),
            "uncertainty": uncertainty,
            "error_classes": list(error_classes),
        },
        "ordinary_Z_order": order,
        "segment": "collar",
        "word": "collar-prior-group",
        "channel": "collar_group",
        "cell": None,
        "receipt": "synthetic://collar",
        "path": group_path,
    }


def _source_index(rows: Mapping[int, Sequence[Mapping[str, Any]]]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for order in ORDERS:
        for record_index, record in enumerate(rows[order]):
            for source_id in record["origin"]["source_ids"]:
                index.setdefault(source_id, {"source_id": source_id, "occurrences": []})["occurrences"].append({
                    "receipt": record["receipt"],
                    "path": record["path"],
                    "ordinary_Z_order": order,
                    "record_index": record_index,
                })
    return index


def _packet(context: Any) -> dict[str, Any]:
    family = {
        "family": "synthetic-n1-segmented-inner",
        "revision": "fixture-v2",
        "owner_identity_proved": False,
    }
    logs = {name: "synthetic:" + name for name in SOURCE_POWER_BASIS}
    shared_source_id = "fixture:shared:core-leading"
    regular_powers = (0, 1, 0, -1)
    uncertain_powers = (1, 0, 0, -1)
    core: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    collar: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    expected_formal: dict[int, dict[tuple[Any, ...], Fraction]] = {0: {}, 1: {}}
    expected_core_group: dict[int, dict[tuple[Any, ...], Fraction]] = {0: {}, 1: {}}
    expected_collar_group: dict[int, dict[tuple[Any, ...], Fraction]] = {0: {}, 1: {}}

    for order in ORDERS:
        if order == 0:
            w0_value, w1_value, core_uncertain, collar_regular, collar_uncertain = (
                Fraction(1, 8), Fraction(1, 16), Fraction(1, 16), Fraction(1, 16), Fraction(-1, 8)
            )
        else:
            w0_value, w1_value, core_uncertain, collar_regular, collar_uncertain = (
                Fraction(-1, 8), Fraction(1, 16), Fraction(-1, 16), Fraction(1, 32), Fraction(1, 8)
            )
        core[order].extend((
            _term(
                context,
                source_id=shared_source_id,
                powers=regular_powers,
                value=w0_value,
                uncertainty=False,
                error_classes=("synthetic:regular",),
                temporal_tag=f"core:W0:Z{order}",
                channel="axial",
                order=order,
                cell=0,
                word="W0",
                receipt="synthetic://core_w0",
                path=["synthetic://core_w0", "records", str(order), "w0"],
            ),
            _term(
                context,
                source_id=shared_source_id,
                powers=regular_powers,
                value=w1_value,
                uncertainty=False,
                error_classes=("synthetic:regular",),
                temporal_tag=f"core:W1-K0:Z{order}",
                channel="axial",
                order=order,
                cell=0,
                word="W1-K0",
                receipt="synthetic://core_w1_k0",
                path=["synthetic://core_w1_k0", "records", str(order), "w1"],
            ),
            _term(
                context,
                source_id=f"fixture:core:uncertain:Z{order}",
                powers=uncertain_powers,
                value=core_uncertain,
                uncertainty=True,
                error_classes=("synthetic:uncertain",),
                temporal_tag=f"core:W0-uncertain:Z{order}",
                channel="swirl",
                order=order,
                cell=1,
                word="W0",
                receipt="synthetic://core_w0",
                path=["synthetic://core_w0", "records", str(order), "uncertain"],
            ),
        ))
        collar[order].extend((
            _prior_group(
                context,
                order=order,
                group_index=0,
                powers=regular_powers,
                group_value=collar_regular,
                member_values=(collar_regular / 2, collar_regular / 2),
                uncertainty=False,
                error_classes=("synthetic:regular",),
            ),
            _prior_group(
                context,
                order=order,
                group_index=1,
                powers=uncertain_powers,
                group_value=collar_uncertain,
                member_values=(collar_uncertain / 2, collar_uncertain / 2),
                uncertainty=True,
                error_classes=("synthetic:uncertain",),
            ),
        ))
        regular_key = (regular_powers, False, ((str, "synthetic:regular"),))
        uncertain_key = (uncertain_powers, True, ((str, "synthetic:uncertain"),))
        core_regular = w0_value + w1_value
        expected_core_group[order] = {regular_key: core_regular, uncertain_key: core_uncertain}
        expected_collar_group[order] = {regular_key: collar_regular, uncertain_key: collar_uncertain}
        expected_formal[order] = {
            regular_key: core_regular + collar_regular,
            uncertain_key: core_uncertain + collar_uncertain,
        }

    collar_source_index = {
        member["record"]["origin"]["source_ids"][0]: {
            "source_id": member["record"]["origin"]["source_ids"][0],
            "path": member["record"]["path"],
        }
        for order in ORDERS
        for group in collar[order]
        for member in group["prior_group_members"]
    }
    return {
        "core": core,
        "collar": collar,
        "family": family,
        "logs": logs,
        "geometry": {
            "rho_seam": True,
            "core_rho_formula": "4*t^2",
            "collar_rho_formula": "4*exp(hb*sc*q/4)",
            "seam_scope": "geometric interface only",
            "endpoint_cancellation_proved": False,
        },
        "shared_source_id": shared_source_id,
        "expected_formal": expected_formal,
        "expected_core_group": expected_core_group,
        "expected_collar_group": expected_collar_group,
        "core_source_index": _source_index(core),
        "collar_source_index": collar_source_index,
    }


def _key(value: Mapping[str, Any]) -> tuple[Any, ...]:
    marker = value.get("group_key", value)
    origin = value.get("origin", {})
    return (
        tuple(marker.get("powers", value.get("powers", ()))),
        bool(marker.get("uncertainty", origin.get("uncertainty"))),
        tuple((str, item) for item in marker.get("error_classes", origin.get("error_classes", ()))),
    )


def _path_key(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple(sorted((str(key), _path_key(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_path_key(item) for item in value)
    return (type(value).__name__, value)


def _member_key(segment: str, order: int, value: Mapping[str, Any]) -> tuple[Any, ...]:
    return (segment, order, str(value.get("receipt")), _path_key(value.get("path")))


def _group_rows(result: Mapping[str, Any], segment: str, order: int) -> Sequence[Any]:
    grouped = result.get("segment_groups", result.get("groups_by_segment", {}))
    if not isinstance(grouped, Mapping) or not isinstance(grouped.get(segment), Mapping):
        return ()
    rows = grouped[segment]
    return rows.get(order, rows.get(str(order), ()))


def _formal_rows(result: Mapping[str, Any], order: int) -> Sequence[Any]:
    formal = result.get("formal_sum", {})
    if not isinstance(formal, Mapping):
        return ()
    return formal.get(order, formal.get(str(order), ()))


def _signed(value: Any) -> tuple[bool, bool]:
    try:
        return bool(value.a > 0), bool(value.b < 0)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False, False


def _validate(context: Any, case: Mapping[str, Any], result: Mapping[str, Any]) -> dict[str, bool]:
    core = case["core"]
    collar = case["collar"]
    expected_records: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    for segment, rows in (("core", core), ("collar", collar)):
        for order in ORDERS:
            for record in rows[order]:
                key = _member_key(segment, order, record)
                if key in expected_records:
                    raise ValueError("synthetic input receipt/path is not unique")
                expected_records[key] = record

    actual_records: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    nested_member_count = 0
    for segment in ("core", "collar"):
        for order in ORDERS:
            for group in _group_rows(result, segment, order):
                if not isinstance(group, Mapping):
                    raise ValueError("candidate group is not a mapping")
                members = group.get("members", ())
                if not isinstance(members, Sequence) or isinstance(members, (str, bytes, bytearray)):
                    raise ValueError("candidate group members are not a sequence")
                for member in members:
                    if not isinstance(member, Mapping):
                        raise ValueError("candidate member is not a mapping")
                    key = _member_key(segment, order, member)
                    if key in actual_records or key not in expected_records:
                        raise ValueError("candidate receipt/path multiplicity mismatch")
                    expected_record = expected_records[key]
                    original = member.get("record")
                    if original is not expected_record:
                        raise ValueError("candidate did not retain the original input record")
                    if member.get("origin") is not expected_record.get("origin"):
                        raise ValueError("candidate changed member origin")
                    if member.get("provenance") is not expected_record.get("origin"):
                        raise ValueError("candidate changed member provenance")
                    if original.get("coefficient") is not expected_record.get("coefficient"):
                        raise ValueError("candidate changed member coefficient")
                    if segment == "collar":
                        prior_members = original.get("prior_group_members", ())
                        if not isinstance(prior_members, list) or len(prior_members) < 2:
                            raise ValueError("prior collar group was not retained as a grouped record")
                        nested_member_count += len(prior_members)
                    actual_records[key] = member
    if set(actual_records) != set(expected_records):
        raise ValueError("candidate omitted or duplicated an input receipt/path")

    result_index = result.get("source_index", {})
    core_index = result_index.get("core", {}) if isinstance(result_index, Mapping) else {}
    collar_index = result_index.get("collar", {}) if isinstance(result_index, Mapping) else {}
    expected_core_ids = set(case["core_source_index"])
    expected_collar_ids = set(case["collar_source_index"])

    formal_checks: list[bool] = []
    for order in ORDERS:
        seen: set[tuple[Any, ...]] = set()
        for term in _formal_rows(result, order):
            if not isinstance(term, Mapping):
                continue
            key = _key(term)
            expected_value = case["expected_formal"][order].get(key)
            formal_checks.append(expected_value is not None)
            if expected_value is None or key in seen:
                continue
            formal_checks.append(_same_interval(term.get("coefficient"), _interval(context, expected_value)))
            seen.add(key)
        formal_checks.append(seen == set(case["expected_formal"][order]))

    core_group_checks: list[bool] = []
    collar_group_checks: list[bool] = []
    collar_direct_checks: list[bool] = []
    for order in ORDERS:
        core_seen: set[tuple[Any, ...]] = set()
        for group in _group_rows(result, "core", order):
            if isinstance(group, Mapping):
                key = _key(group)
                expected_value = case["expected_core_group"][order].get(key)
                core_group_checks.append(expected_value is not None)
                if expected_value is not None:
                    core_group_checks.append(_same_interval(group.get("coefficient"), _interval(context, expected_value)))
                    core_seen.add(key)
        core_group_checks.append(core_seen == set(case["expected_core_group"][order]))
        collar_seen: set[tuple[Any, ...]] = set()
        for group in _group_rows(result, "collar", order):
            if isinstance(group, Mapping):
                key = _key(group)
                expected_value = case["expected_collar_group"][order].get(key)
                collar_group_checks.append(expected_value is not None)
                if expected_value is not None:
                    collar_group_checks.append(_same_interval(group.get("coefficient"), _interval(context, expected_value)))
                    collar_seen.add(key)
                for member in group.get("members", ()):
                    if isinstance(member, Mapping) and isinstance(member.get("record"), Mapping):
                        collar_direct_checks.append(member["record"] is next(
                            (record for record in collar[order] if _member_key("collar", order, record) == _member_key("collar", order, member)),
                            None,
                        ))
        collar_group_checks.append(collar_seen == set(case["expected_collar_group"][order]))

    output_coefficients = [
        term.get("coefficient")
        for order in ORDERS
        for term in _formal_rows(result, order)
        if isinstance(term, Mapping)
    ]
    positive = any(_signed(value)[0] for value in output_coefficients)
    negative = any(_signed(value)[1] for value in output_coefficients)
    geometry = result.get("geometry_certificate")
    scope = result.get("scope", {})
    return {
        "source_family_exact": result.get("source_family") == case["family"],
        "source_logs_exact": result.get("source_logs") == case["logs"],
        "source_atoms_exact": result.get("source_atoms") == case["logs"],
        "C0_present": bool(_formal_rows(result, 0)),
        "Z1_present": bool(_formal_rows(result, 1)),
        "signed_nonzero_formal_output": positive and negative,
        "same_source_id_in_W0_and_W1_allowed": sum(
            1 for order in ORDERS for record in core[order]
            if case["shared_source_id"] in record["origin"]["source_ids"]
        ) == 4,
        "receipt_path_unique_and_exactly_once": set(actual_records) == set(expected_records),
        "core_members_retained_once": sum(
            len(group.get("members", ()))
            for order in ORDERS
            for group in _group_rows(result, "core", order)
            if isinstance(group, Mapping)
        ) == sum(len(core[order]) for order in ORDERS),
        "collar_prior_groups_not_reexpanded": (
            sum(len(collar[order]) for order in ORDERS) == 4
            and sum(
                len(group.get("members", ()))
                for order in ORDERS
                for group in _group_rows(result, "collar", order)
                if isinstance(group, Mapping)
            ) == 4
        ),
        "prior_nested_members_retained": nested_member_count == 8,
        "core_source_index_complete": isinstance(core_index, Mapping) and expected_core_ids.issubset(set(core_index)),
        "collar_source_index_complete": isinstance(collar_index, Mapping) and expected_collar_ids.issubset(set(collar_index)),
        "different_error_key_groups_retained": all(
            len(_group_rows(result, segment, order)) == 2
            for segment in ("core", "collar")
            for order in ORDERS
        ),
        "core_group_coefficients_exact": all(core_group_checks),
        "collar_group_coefficients_exact": all(collar_group_checks),
        "prior_group_coefficients_consumed_directly": all(collar_direct_checks) and len(collar_direct_checks) == 4,
        "formal_core_plus_collar_coefficients_exact": all(formal_checks),
        "rho_seam_is_geometric_only": (
            isinstance(geometry, Mapping)
            and geometry.get("rho_seam") is True
            and geometry.get("endpoint_cancellation_proved") is False
        ),
        "w1_scope_is_partial": isinstance(scope, Mapping) and scope.get("core", {}).get("full_original_core_W1_missing") is True,
        "P1_is_not_moment": result.get("P1_is_moment") is False and result.get("P1_mixed") in (False, None),
        "physical_scale_not_reapplied": result.get("physical_scale_reapplied") is False and result.get("physical_factors_applied") in (False, None),
        "no_closure_or_endpoint_claim": result.get("endpoint_cancellation_proved") is False and result.get("terminal_closure_proved") is False,
        "no_ledger_change": result.get("ledger_212_changed") is False and result.get("ledger_appended") is False and result.get("writes_ledger") is False,
        "owner_identity_not_proved": result.get("owner_identity_proved") is False,
        "same_interval_context_left_to_caller": result.get("same_interval_context_validated") is False,
    }


def run_fixture(leaf: Any, context: Any) -> dict[str, Any]:
    """Call the supplied candidate exactly once on the synthetic packet."""

    case = _packet(context)
    checks: dict[str, bool] = {
        "leaf_build_available": callable(getattr(leaf, "build", None)),
        "exact_rational_binary_fixture": True,
        "collar_has_two_members_per_prior_group": all(
            len(group["prior_group_members"]) >= 2
            for order in ORDERS
            for group in case["collar"][order]
        ),
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
            "implementation_imported": False,
            "leaf_calls": 1,
        }
    try:
        result = build(
            context,
            core=case["core"],
            collar=case["collar"],
            source_family=case["family"],
            source_atoms=case["logs"],
            source_logs=case["logs"],
            geometry_certificate=case["geometry"],
            owner_certificate={"source_family_exact": True, "rho_seam_only": True, "owner_identity_proved": False},
            scale_certificate={"already_physically_scaled": True, "physical_scale_reapplied": False},
        )
        checks["build_completed"] = isinstance(result, Mapping)
        if isinstance(result, Mapping):
            checks.update(_validate(context, case, result))
    except BaseException as exc:
        checks["build_completed"] = False
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "check_count": len(checks),
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
            "leaf_calls": 1,
        }
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "exact_rational_inputs": True,
        "signed_nonzero_inputs": True,
        "prior_group_count": 4,
        "prior_group_member_count": 8,
        "prior_groups_not_reexpanded": True,
        "implementation_imported": False,
        "old_code_called": False,
        "leaf_calls": 1,
    }


def run(leaf: Any, context: Any) -> dict[str, Any]:
    return run_fixture(leaf, context)
