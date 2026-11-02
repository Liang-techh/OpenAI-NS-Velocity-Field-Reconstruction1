"""Independent v4 synthetic audit for :mod:`w0_joined_jets_v4`.

The fixture uses a new source family, IDs, paths and rational intervals.  It
checks the direct product-state shape as well as source-bound rejection paths.
It performs no source evaluation, quadrature, interval summation or receipt
access.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from typing import Any


FIELDS = ("F1", "Uz1", "K1")
OLD_FIELDS = FIELDS + ("d_x_F1",)
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1, 2)
CELL_COUNT = 4
COMPONENT_IDS = {"F1": 1, "Uz1": 2, "K1": 3, "P1": 4, "d_x_F1": 5, "d_x_Uz1": 6}
ROWS = {"F1": "5", "Uz1": "6", "K1": "6", "d_x_F1": "5"}
SOURCE_FAMILY = "fixture-v4:joined-w0:family-7"
SOURCE_ATOMS = {name: f"fixture-v4:atom:{index + 11}" for index, name in enumerate(ATOMS)}
ALLOWED_STATE = {
    (component, order)
    for order in (0, 1)
    for component in (1, 2, 3, 5)
} | {(1, 2), (2, 2), (3, 2)}


def _point(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _interval(c: Any, lower: Fraction, upper: Fraction) -> Any:
    left = _point(c, lower)
    right = _point(c, upper)
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _cells(c: Any) -> list[dict[str, Any]]:
    sc = _point(c, Fraction(7, 16384))
    return [
        {
            "index": index,
            "q_lower": _point(c, Fraction(index, 4)),
            "q_upper": _point(c, Fraction(index + 1, 4)),
            "s_lower": sc * index / 16,
            "s_upper": sc * (index + 1) / 16,
        }
        for index in range(CELL_COUNT)
    ]


def _term(
    c: Any,
    *,
    role: str,
    order: int,
    field: str,
    source_id: str,
    cell: int | None = None,
    lower: Fraction,
    upper: Fraction,
) -> dict[str, Any]:
    path = f"fixture-v4/{role}/record/{source_id}/Z{order}/{field}"
    errors = ["fixture-v4:error-bound", f"fixture-v4:{role}:error"]
    ref: dict[str, Any] = {
        "source_id": source_id,
        "source_path": path,
        "packet_role": role,
    }
    origin: dict[str, Any] = {
        "source_ids": [source_id],
        "source_paths": [path],
        "source_node_refs": [ref],
        "error_classes": errors,
        "uncertainty": order == 1,
        "packet_role": role,
        "ordinary_Z_order": order,
        "field": field,
        "source_row": ROWS[field],
        "operation": "fixture_v4_join_source",
    }
    if cell is not None:
        origin["phase_cell"] = cell
        ref["phase_cell"] = cell
    term = {
        "powers": [2 + order, (order + len(field) + 1) % 4, -1 - order, 2],
        "coefficient": _interval(c, lower, upper),
        "source_ids": [source_id],
        "source_paths": [path],
        "source_node_refs": [ref],
        "uncertainty": order == 1,
        "error_classes": errors,
        "origin": origin,
    }
    if cell is not None:
        term["phase_cell"] = cell
    if order == 2:
        term["source_z_derivative_is_ordinary"] = True
        term["ordinary_Z_factorial_divided"] = False
    return term


def _payload(c: Any, role: str) -> dict[str, Any]:
    orders = (0, 1) if role == "z0_z1" else (2,)
    fields = OLD_FIELDS if role == "z0_z1" else FIELDS
    core: dict[Any, Any] = {}
    for order in orders:
        row: dict[str, Any] = {}
        for index, field in enumerate(fields):
            base = Fraction((order + 3) * (index + 2), 53)
            row[field] = [_term(
                c,
                role=role,
                order=order,
                field=field,
                source_id=f"v4-{role}-core-{order}-{field}-q9",
                lower=-base - Fraction(1, 211) if index % 2 else base,
                upper=base + Fraction(2, 211),
            )]
        core[order] = row
    if role == "z2":
        core = core[2]

    if role == "z0_z1":
        local: dict[int, Any] = {}
        for order in orders:
            local[order] = {}
            for cell in range(CELL_COUNT):
                local[order][cell] = {}
                for index, field in enumerate(fields):
                    base = Fraction((cell + 2) * (index + order + 2), 97)
                    local[order][cell][field] = [_term(
                        c,
                        role=role,
                        order=order,
                        field=field,
                        cell=cell,
                        source_id=f"v4-{role}-local-{order}-{cell}-{field}-q9",
                        lower=-base if (cell + index) % 2 else base,
                        upper=base + Fraction(1, 239),
                    )]
    else:
        local = {}
        for cell in range(CELL_COUNT):
            local[cell] = {}
            for index, field in enumerate(fields):
                base = Fraction((cell + 4) * (index + 1), 113)
                local[cell][field] = [_term(
                    c,
                    role=role,
                    order=2,
                    field=field,
                    cell=cell,
                    source_id=f"v4-{role}-local-2-{cell}-{field}-q9",
                    lower=-base - Fraction(1, 263) if index == 2 else base,
                    upper=base + Fraction(3, 263),
                )]
    return {
        "schema": "fixture-v4:source-packet",
        "source_family": SOURCE_FAMILY,
        "source_atoms": deepcopy(SOURCE_ATOMS),
        "source_power_basis": list(ATOMS),
        "ordinary_Z_orders": list(orders),
        "source_z2_is_ordinary_derivative": role == "z2",
        "ordinary_Z_factorial_divided": False,
        "fields": list(fields),
        "core_transport": core,
        "local": local,
        "q_cells": _cells(c),
    }


def _packet(c: Any, role: str) -> dict[str, Any]:
    return {
        "source_verified": True,
        "verified_source_package": True,
        "packet_role": role,
        "source_family": SOURCE_FAMILY,
        "source_atoms": deepcopy(SOURCE_ATOMS),
        "source_power_basis": list(ATOMS),
        "payload": _payload(c, role),
    }


def _same_term(left: Any, right: Any) -> bool:
    return (
        left.get("powers") == right.get("powers")
        and _same_interval(left.get("coefficient"), right.get("coefficient"))
        and left.get("source_ids") == right.get("source_ids")
        and left.get("source_paths") == right.get("source_paths")
        and left.get("source_node_refs") == right.get("source_node_refs")
        and left.get("uncertainty") == right.get("uncertainty")
        and left.get("error_classes") == right.get("error_classes")
        and left.get("source_z_derivative_is_ordinary") == right.get("source_z_derivative_is_ordinary")
        and left.get("ordinary_Z_factorial_divided") == right.get("ordinary_Z_factorial_divided")
        and left.get("origin") == right.get("origin")
    )


def _reject(build: Any, old: Any, new: Any) -> bool:
    try:
        build(old, new)
    except (KeyError, TypeError, ValueError):
        return True
    return False


def run(leaf: Any, c: Any) -> dict[str, Any]:
    build = leaf if callable(leaf) else leaf.build
    old = _packet(c, "z0_z1")
    new = _packet(c, "z2")
    try:
        result = build(old, new)
    except (KeyError, TypeError, ValueError) as exc:
        return {
            "schema": "w0_joined_jets_fixture_v4",
            "passed": False,
            "check_count": 1,
            "checks": {"valid_join": False},
            "synthetic_only": True,
            "error": f"{type(exc).__name__}: {exc}",
        }

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "w0_joined_jets_v4",
        "source_family": result.get("source_family") == SOURCE_FAMILY,
        "source_bound": result.get("source_bound") is True,
        "basis": tuple(result.get("source_power_basis", ())) == ATOMS,
        "three_orders": tuple(result.get("ordinary_Z_orders", ())) == ORDERS,
        "four_cells": isinstance(result.get("q_cells"), list) and len(result["q_cells"]) == CELL_COUNT,
        "core_history": result.get("core_history_reset") is False and result.get("core_retained_once") is True,
        "separate_sectors": result.get("local_terms_are_separate") is True,
        "source_flags": all(result.get(name) is True for name in (
            "source_ids_preserved", "intervals_preserved", "formal_powers_preserved", "error_classes_preserved"
        )),
        "no_evaluation": all(result.get(name) is False for name in (
            "source_values_computed", "actual_forced_evaluation", "quadrature_used", "formal_powers_materialized"
        )),
        "z2_ordinary_no_factorial": result.get("source_z2_is_ordinary_derivative") is True
        and result.get("ordinary_Z_factorial_divided") is False,
        "no_g5_or_ledger": result.get("g5_added") is False and result.get("g5_recomputed") is False
        and result.get("m15_computed") is False and result.get("m15_ledger_changed") is False
        and result.get("m15_terms_added") == 0 and result.get("m15_ledger_size_before") == result.get("m15_ledger_size_after") == 212,
        "feedback_retained": result.get("feedback_remainder", {}).get("eliminated") is False
        and result.get("feedback_remainder", {}).get("zero_padded") is False,
    }

    core = result.get("core_transport")
    local = result.get("local")
    state = result.get("state")
    unknown = result.get("unknown_state")
    checks["canonical_core_shape"] = isinstance(core, dict) and set(core) == {0, 1, 2}
    checks["canonical_local_shape"] = isinstance(local, dict) and set(local) == {0, 1, 2}
    checks["core_count_not_repeated"] = checks["canonical_core_shape"] and (
        sum(len(core[o][field]) for o in (0, 1) for field in OLD_FIELDS)
        + sum(len(core[2][field]) for field in FIELDS)
        == 11
    )
    checks["core_records_exact"] = checks["canonical_core_shape"] and all(
        _same_term(
            core[o][field][0],
            _payload(c, "z0_z1")["core_transport"][o][field][0],
        )
        for o in (0, 1)
        for field in OLD_FIELDS
    ) and all(
        _same_term(core[2][field][0], _payload(c, "z2")["core_transport"][field][0])
        for field in FIELDS
    )
    checks["local_provenance_exact"] = checks["canonical_local_shape"] and all(
        local[order][cell][field][0].get("origin", {}).get("phase_cell") == cell
        and str(local[order][cell][field][0].get("origin", {}).get("source_row")) == ROWS[field]
        for order in ORDERS
        for cell in range(CELL_COUNT)
        for field in (OLD_FIELDS if order < 2 else FIELDS)
    )
    checks["local_error_metadata"] = checks["canonical_local_shape"] and all(
        term.get("error_classes") and term.get("origin", {}).get("error_classes")
        for order in ORDERS
        for cell in range(CELL_COUNT)
        for field in (OLD_FIELDS if order < 2 else FIELDS)
        for term in local[order][cell][field]
    )

    checks["product_state_shape"] = isinstance(state, dict) and set(state) == {0, 1, 2, 3}
    if checks["product_state_shape"]:
        checks["product_keys_allowed"] = all(set(state[cell]) == ALLOWED_STATE for cell in state)
        checks["product_terms_are_records"] = all(
            isinstance(state[cell][key][sector], list)
            and state[cell][key][sector]
            and isinstance(state[cell][key][sector][0], dict)
            and "coefficient" in state[cell][key][sector][0]
            and state[cell][key][sector][0].get("kind") != "source_term_ref"
            for cell in state
            for key in state[cell]
            for sector in ("core", "local")
        )
        checks["product_terms_canonical"] = all(
            state[cell][key]["core"] is core[key[1]][{1: "F1", 2: "Uz1", 3: "K1", 5: "d_x_F1"}[key[0]]]
            and state[cell][key]["local"] is local[key[1]][cell][{1: "F1", 2: "Uz1", 3: "K1", 5: "d_x_F1"}[key[0]]]
            for cell in state
            for key in state[cell]
        )
        checks["no_implicit_sum"] = all(
            state[cell][key].get("implicit_addition") is False
            and state[cell][key].get("core_refs")
            and state[cell][key].get("local_refs")
            for cell in state
            for key in state[cell]
        )
    else:
        checks["product_keys_allowed"] = False
        checks["product_terms_are_records"] = False
        checks["product_terms_canonical"] = False
        checks["no_implicit_sum"] = False

    checks["unknowns_explicit"] = isinstance(unknown, dict) and all(
        unknown[cell][(component, order)]["status"] == "unknown"
        for cell in range(CELL_COUNT)
        for order in ORDERS
        for component in (4, 6)
    ) and all(unknown[cell][(5, 2)]["status"] == "unknown" for cell in range(CELL_COUNT))

    missing = deepcopy(old)
    del missing["payload"]["local"][1][2]["K1"]
    checks["missing_data_rejected"] = _reject(build, missing, new)
    fake_owner = deepcopy(new)
    fake_owner["payload"]["local"][0]["F1"][0]["packet_role"] = "z0_z1"
    checks["pseudo_ownership_rejected"] = _reject(build, old, fake_owner)
    wrong_cell = deepcopy(new)
    wrong_cell["payload"]["q_cells"][3]["q_upper"] = _point(c, Fraction(1, 2))
    checks["wrong_q_cell_rejected"] = _reject(build, old, wrong_cell)
    unverified = deepcopy(new)
    unverified["source_verified"] = False
    unverified["verified_source_package"] = False
    checks["unverified_rejected"] = _reject(build, old, unverified)
    family_mismatch = deepcopy(new)
    family_mismatch["source_family"] = "fixture-v4:other-family"
    checks["family_mismatch_rejected"] = _reject(build, old, family_mismatch)
    wrong_atoms = deepcopy(new)
    wrong_atoms["payload"]["source_atoms"]["hb"] = "fixture-v4:wrong-atom"
    checks["wrapper_payload_atoms_mismatch_rejected"] = _reject(build, old, wrong_atoms)
    wrong_factorial = deepcopy(new)
    wrong_factorial["payload"]["ordinary_Z_factorial_divided"] = True
    checks["factorial_divided_rejected"] = _reject(build, old, wrong_factorial)

    return {
        "schema": "w0_joined_jets_fixture_v4",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "new_synthetic_family": SOURCE_FAMILY,
        "actual_source_values_computed": False,
        "quadrature_used": False,
        "ordinary_orders": list(ORDERS),
        "cell_count": CELL_COUNT,
    }


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    return run(build, c)


__all__ = ["ATOMS", "FIELDS", "OLD_FIELDS", "ROWS", "run", "run_fixture"]
