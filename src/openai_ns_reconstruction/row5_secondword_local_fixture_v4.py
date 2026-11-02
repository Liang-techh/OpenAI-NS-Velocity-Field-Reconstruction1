"""Fresh synthetic checks for the factorized local second-word leaf.

The fixture passes a saved-ledger-shaped mapping directly to the leaf; it does
not import or call the v4 producer.  Ordinary inputs are exact directed MPF
endpoint tuples.  Checks cover signed positive/negative K0 prefixes, scalar
application exactly once, source-family isolation, the unknown core inlet,
formal-hb preservation, cell-independent grouping, and scale factors.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "row5_secondword_local_fixture_v4"
INPUT_SCHEMA = "row5_firstword_product_v4"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
SOURCE_ROWS = ("h5", "h5_Z1")
OUTPUT_FIELDS = ("F1", "d_x_F1")

__all__ = ["SCHEMA", "run_fixture", "run"]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None and len(raw) == 2:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpi_", None)
        if nested is not None and len(nested) == 2:
            return nested[1 if upper else 0]
        nested = getattr(endpoint, "_mpf_", None)
        if nested is not None:
            return nested
        return endpoint
    return value


def _raw_endpoints(value: Any) -> tuple[Any, Any]:
    return _raw_endpoint(value, False), _raw_endpoint(value, True)


def _q(c: Any, value: Fraction | int) -> Any:
    rational = Fraction(value)
    return c.mpf(rational.numerator) / rational.denominator


def _payload(c: Any, lower: Fraction | int, upper: Fraction | int) -> dict[str, Any]:
    lower_raw = _raw_endpoint(_q(c, lower), False)
    upper_raw = _raw_endpoint(_q(c, upper), True)
    return {
        "lower_exact_mpf_tuple": tuple(lower_raw),
        "upper_exact_mpf_tuple": tuple(upper_raw),
        "directed": True,
        "function_value": False,
    }


def _decode(c: Any, payload: Mapping[str, Any]) -> Any:
    return c.make_mpf((tuple(payload["lower_exact_mpf_tuple"]), tuple(payload["upper_exact_mpf_tuple"])))


def _same_interval(c: Any, left: Any, right: Any) -> bool:
    left_iv = _decode(c, left) if isinstance(left, Mapping) else left
    right_iv = _decode(c, right) if isinstance(right, Mapping) else right
    return _raw_endpoints(left_iv) == _raw_endpoints(right_iv)


def _node(
    c: Any,
    node_id: str,
    family: str,
    coefficient: Fraction | int,
    powers: Sequence[int] = (0, 0, 0, 0),
    *,
    scalar: int = 1,
    error: bool = False,
) -> dict[str, Any]:
    return {
        "node_id": node_id,
        "source_family": family,
        "coefficient_enclosure": _payload(c, coefficient, coefficient),
        "formal_powers": list(powers),
        "source_ids": [node_id],
        "source_paths": [["fixture_v4", "source_nodes", node_id]],
        "source_error_classes": ["fixture_v4:error:" + node_id] if error else [],
        "error_classes": ["fixture_v4:error:" + node_id] if error else [],
        "error_only": error,
        "uncertainty": error,
        "symmetric_error": error,
        "scalar": scalar,
    }


def _product(
    c: Any,
    product_id: str,
    row: str,
    source_node_id: str,
    coefficient: Fraction | int,
    powers: Sequence[int] = (0, 0, 0, 0),
    *,
    scalar: int = 1,
    state_node_id: str | None = None,
    error: bool = False,
) -> dict[str, Any]:
    return {
        "product_id": product_id,
        "row": row,
        "source_node_id": source_node_id,
        "state_node_id": state_node_id,
        "source_only": False,
        "scalar": scalar,
        "formal_powers": list(powers),
        "coefficient_enclosure": _payload(c, coefficient, coefficient),
        "scalar_applied_to_coefficient": False,
        "term_semantics": "scalar * coefficient_enclosure * formal_monomial; scalar deferred",
        "error_only": error,
        "uncertainty": error,
        "symmetric_error": error,
    }


def _packet(c: Any) -> dict[str, Any]:
    source_nodes = {
        "v4node:positive": _node(c, "v4node:positive", "C0", Fraction(5, 3), (0, 2, 0, 0)),
        "v4node:negative": _node(c, "v4node:negative", "B1", Fraction(-7, 5), (1, 1, 0, 0)),
        "v4node:scalar": _node(c, "v4node:scalar", "C0", Fraction(11, 4), (0, 0, 2, 0), scalar=2),
        "v4node:error": _node(c, "v4node:error", "B1", -3, (0, 0, 0, 2)),
        "v4node:huge": _node(c, "v4node:huge", "C0", 10**7, (8, -4, 0, 0)),
        "v4node:g5": _node(c, "v4node:g5", "g5", 47, error=True),
        "v4node:axis": _node(c, "v4node:axis", "axis", 43, scalar=2, error=True),
        "v4node:z1": _node(c, "v4node:z1", "C0", Fraction(-5, 2)),
        "v4node:z1_error": _node(c, "v4node:z1_error", "B1", Fraction(13, 6), error=True),
        "v4node:z1_g5": _node(c, "v4node:z1_g5", "g5", 31),
        "v4node:z1_axis": _node(c, "v4node:z1_axis", "axis", 37, scalar=2),
    }
    rows = {
        "h5": [
            _product(c, "v4p:positive", "h5", "v4node:positive", Fraction(5, 3), (0, 2, 0, 0), state_node_id="state:V4:W1:Z0:0"),
            _product(c, "v4p:negative", "h5", "v4node:negative", Fraction(-7, 5), (1, 1, 0, 0), state_node_id="state:V4:W1:Z0:0"),
            _product(c, "v4p:scalar", "h5", "v4node:scalar", Fraction(11, 4), (0, 0, 2, 0), scalar=2, state_node_id="state:V4:W1:Z0:0"),
            _product(c, "v4p:error", "h5", "v4node:error", -3, (0, 0, 0, 2), state_node_id="state:V4:W1:Z0:error"),
            _product(c, "v4p:huge", "h5", "v4node:huge", 10**7, (8, -4, 0, 0), state_node_id="state:V4:W1:Z0:0"),
            _product(c, "v4p:g5", "h5", "v4node:g5", 47, error=True),
            _product(c, "v4p:axis", "h5", "v4node:axis", 43, scalar=2, error=True),
        ],
        "h5_Z1": [
            _product(c, "v4p:z1", "h5_Z1", "v4node:z1", Fraction(-5, 2), state_node_id="state:V4:W1:Z0:0"),
            _product(c, "v4p:z1_error", "h5_Z1", "v4node:z1_error", Fraction(13, 6), state_node_id="state:V4:W1:Z0:error", error=True),
            _product(c, "v4p:z1_g5", "h5_Z1", "v4node:z1_g5", 31),
            _product(c, "v4p:z1_axis", "h5_Z1", "v4node:z1_axis", 37, scalar=2),
        ],
    }
    return {
        "schema": INPUT_SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "source_nodes": source_nodes,
        "state_nodes": {
            "state:V4:W1:Z0:0": {
                "source_ids": ["state:V4:W1:Z0:0"],
                "source_paths": [["fixture_v4", "state", "W1", "Z0", 0]],
                "coefficient_enclosure": _payload(c, 1, 1),
                "formal_powers": [0, 0, 0, 0],
                "error_only": False,
                "uncertainty": False,
                "symmetric_error": False,
            },
            "state:V4:W1:Z0:error": {
                "source_ids": ["state:V4:W1:Z0:error"],
                "source_paths": [["fixture_v4", "state", "W1", "Z0", "error"]],
                "coefficient_enclosure": _payload(c, 1, 1),
                "formal_powers": [0, 0, 0, 0],
                "error_only": True,
                "uncertainty": True,
                "symmetric_error": True,
            }
        },
        "rows": rows,
        "formal_atoms_materialized": False,
        "physical_h5_materialized": False,
    }


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError, ZeroDivisionError):
        return True
    return False


def _callable(value: Any, name: str) -> Any:
    candidate = getattr(value, name, value)
    if not callable(candidate):
        raise TypeError(f"{name} must be callable or expose {name}")
    return candidate


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    build_fn = _callable(build, "build")
    evaluate_fn = getattr(build_fn, "evaluate_term", None)
    sc = _payload(c, Fraction(1, 96), Fraction(1, 96))
    q_cells = [
        _payload(c, 0, Fraction(1, 32)),
        _payload(c, Fraction(1, 32), Fraction(1, 4)),
        _payload(c, Fraction(1, 4), Fraction(3, 4)),
        _payload(c, Fraction(3, 4), 1),
    ]
    packet = _packet(c)
    try:
        result = build_fn(c, packet, hb_admissibility=True, sc=sc, q_cells=q_cells)
    except Exception as exc:
        import traceback
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks: dict[str, bool] = {}
    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "row5_secondword_local_v4"
    checks["input_schema"] = result.get("input_schema") == INPUT_SCHEMA
    checks["four_cells"] = len(result.get("cell_factors", {})) == 4 and result.get("geometry", {}).get("q_cell_count") == 4
    checks["v4_custom_q_geometry"] = _same_interval(c, result["cell_factors"]["q0"]["q_interval"], q_cells[0]) and _same_interval(c, result["cell_factors"]["q1"]["q_interval"], q_cells[1])
    checks["source_rows"] = result.get("source_rows") == {"Z0": "h5", "Z1": "h5_Z1"}
    checks["whole_collar_source_scope"] = (
        result.get("scope", {}).get("source_rows_are_whole_collar_enclosures") is True
        and result.get("scope", {}).get("whole_collar_source_enclosure_verified") is False
    )
    checks["fields"] = all(field in result.get("fields", {}) for field in OUTPUT_FIELDS) and all(
        set(result["fields"][field]) == {"Z0", "Z1"}
        and all(
            result["fields"][field][f"Z{order}"].get("factorized") is True
            and result["fields"][field][f"Z{order}"].get("term_expansion_performed") is False
            and len(result["fields"][field][f"Z{order}"].get("cell_factor_ids", [])) == 4
            for order in (0, 1)
        )
        for field in OUTPUT_FIELDS
    )
    checks["compact_factorized"] = all(
        "product_group_ids" in result["fields"][field][f"Z{order}"]
        and "factor_refs" in result["fields"][field][f"Z{order}"]
        and "terms" not in result["fields"][field][f"Z{order}"]
        for field in OUTPUT_FIELDS
        for order in (0, 1)
    ) and all(
        "cell_factor_id" in bound and "member_product_group_id" in bound
        for bound in result["absolute_coefficient_bounds"]
    )
    checks["product_table_once"] = set(result["products"]) == {
        "v4p:positive",
        "v4p:negative",
        "v4p:scalar",
        "v4p:error",
        "v4p:huge",
        "v4p:z1",
        "v4p:z1_error",
    }
    checks["source_refs_preserved"] = (
        result["products"]["v4p:positive"]["source_node_ref"] == "v4node:positive"
        and result["products"]["v4p:positive"]["state_node_ref"] == "state:V4:W1:Z0:0"
        and _same_interval(c, result["products"]["v4p:positive"]["original_coefficient_enclosure"], packet["rows"]["h5"][0]["coefficient_enclosure"])
        and "v4node:g5" in result["source_nodes"]
        and "v4node:axis" in result["source_nodes"]
    )
    excluded_ids = {item["product_id"] for item in result["excluded_products"]}
    checks["g5_axis_isolation"] = excluded_ids == {"v4p:g5", "v4p:axis", "v4p:z1_g5", "v4p:z1_axis"} and not (
        excluded_ids & set(result["products"])
    )
    checks["not_algebraic_cancellation"] = all(
        item.get("algebraic_cancellation") is False for item in result["excluded_products"]
    ) and result["scope"].get("axis_skip_is_algebraic_cancellation") is False
    checks["core_unknown"] = result["core_feedback_inlet"] == {
        "status": "unknown",
        "known": False,
        "included": False,
        "reset_to_zero": False,
        "scope": "local feedback source only",
    }
    checks["formal_hb_kept"] = (
        result["scope"].get("formal_hb_materialized") is False
        and result["geometry"].get("hb_not_materialized") is True
        and all("hb_value" not in cell for cell in result["cell_factors"].values())
        and result["cell_factors"]["q0"]["factors"]["F1"]["formal_hb_power"] == 2
        and result["cell_factors"]["q0"]["factors"]["d_x_F1"]["formal_hb_power"] == 1
    )
    checks["source_identity_gates"] = all(
        result["source_identity_conditions"].get(key) == "caller condition; not proved here"
        for key in (
            "same_source_uniform_on_each_full_collar_cell",
            "common_source_owner",
            "fixed_normalized_x_geometry",
            "whole_collar_callback_identity",
        )
    )

    product_groups = {group["product_group_id"]: group for group in result["product_groups"]}
    checks["product_groups_exact"] = all(
        group.get("member_count") == len(group.get("member_product_ids", []))
        and all(product_id in result["products"] for product_id in group.get("member_product_ids", []))
        for group in result["product_groups"]
    )
    checks["cell_groups_independent"] = all(
        bound.get("cell_factor_id") in result["cell_factors"]
        and bound.get("member_product_group_id") in product_groups
        and bound.get("member_product_ids_ref", "").endswith(".member_product_ids")
        and bound.get("same_source_reused_across_cells") is True
        and bound.get("cell_bound_is_not_independent_source_cell_data") is True
        for bound in result["absolute_coefficient_bounds"]
    )

    # Independent closed-form kernel masses on one nontrivial synthetic
    # target, including the formal hb factor.  These scalar sample checks
    # complement, but do not replace, the directed source enclosures.
    import mpmath as mp
    mp.mp.dps = 100
    sample_hb = mp.mpf(2) / 5
    sample_sc = mp.mpf(1) / 96
    sample_q = mp.mpf(1) / 32
    sample_s = sample_sc * sample_q / 4
    sample_x = 2 * mp.exp(sample_hb * sample_s / 2)
    sample_S = (sample_sc / 4) * sample_q * mp.exp(sample_sc / 16)
    direct_F = sample_x**2 / 8 - 1 + 2 / sample_x**2
    direct_Fx = sample_x / 4 - 4 / sample_x**3
    checks["closed_form_F_mass_within_bound"] = bool(0 < direct_F <= sample_hb**2 * sample_S**2 / 2)
    checks["closed_form_Fx_mass_within_bound"] = bool(0 < direct_Fx <= sample_hb * sample_S)

    sc_iv = _decode(c, sc)
    qhi = _q(c, Fraction(1, 32))
    S_i = (sc_iv / 4) * qhi * c.exp(sc_iv / 16)
    K0_prefix_coefficient = c.make_mpf((_raw_endpoint(_q(c, 0), False), _raw_endpoint(S_i * S_i / 2, True)))
    d_x_prefix_coefficient = c.make_mpf((_raw_endpoint(_q(c, 0), False), _raw_endpoint(S_i, True)))
    expected_positive = _q(c, Fraction(5, 3)) * K0_prefix_coefficient
    expected_negative = _q(c, Fraction(-7, 5)) * K0_prefix_coefficient
    expected_scalar = _q(c, 2 * Fraction(11, 4)) * K0_prefix_coefficient
    expected_scalar_dx = _q(c, 2 * Fraction(11, 4)) * d_x_prefix_coefficient
    if not callable(evaluate_fn):
        checks["evaluate_term_exported"] = False
        positive_eval = negative_eval = scalar_eval = scalar_dx_eval = {}
    else:
        checks["evaluate_term_exported"] = True
        positive_eval = evaluate_fn(c, result["products"]["v4p:positive"], result["cell_factors"]["q0"]["factors"]["F1"])
        negative_eval = evaluate_fn(c, result["products"]["v4p:negative"], result["cell_factors"]["q0"]["factors"]["F1"])
        scalar_eval = evaluate_fn(c, result["products"]["v4p:scalar"], result["cell_factors"]["q0"]["factors"]["F1"])
        scalar_dx_eval = evaluate_fn(c, result["products"]["v4p:scalar"], result["cell_factors"]["q0"]["factors"]["d_x_F1"])
    checks["positive_direct_K0_prefix"] = callable(evaluate_fn) and _same_interval(c, positive_eval["signed_coefficient_enclosure"], expected_positive)
    checks["negative_direct_K0_prefix"] = callable(evaluate_fn) and _same_interval(c, negative_eval["signed_coefficient_enclosure"], expected_negative)
    checks["scalar_once_F1"] = callable(evaluate_fn) and _same_interval(c, scalar_eval["signed_coefficient_enclosure"], expected_scalar)
    checks["scalar_once_dx"] = callable(evaluate_fn) and _same_interval(c, scalar_dx_eval["signed_coefficient_enclosure"], expected_scalar_dx)
    checks["ordinary_Z1_source"] = result["fields"]["F1"]["Z1"]["source_row"] == "h5_Z1"
    checks["negative_sign_retained"] = callable(evaluate_fn) and _decode(c, negative_eval["signed_coefficient_enclosure"]).a < c.mpf(0)
    checks["evaluate_term_metadata"] = callable(evaluate_fn) and scalar_eval["scalar"] == 2 and scalar_eval["scalar_applied_once"] is True and scalar_eval["formal_hb_not_materialized"] is True
    checks["error_separate"] = (
        any(
            group["error_branch"] == "error"
            and "v4p:error" in product_groups[group["member_product_group_id"]]["member_product_ids"]
            for group in result["absolute_coefficient_bounds"]
        )
        and all(
            "v4p:error" not in product_groups[group["member_product_group_id"]]["member_product_ids"]
            for group in result["absolute_coefficient_bounds"]
            if group["error_branch"] == "finite"
        )
    )
    checks["no_cancellation"] = all(
        group.get("no_cancellation") is True for group in result["absolute_coefficient_bounds"]
    ) and all(group.get("no_cancellation") is True for group in result["product_groups"])
    huge_group = next(
        (
            group
            for group in result["absolute_coefficient_bounds"]
            if "v4p:huge" in product_groups[group["member_product_group_id"]]["member_product_ids"]
            and group["field"] == "F1"
            and group["ordinary_Z_order"] == 0
        ),
        None,
    )
    checks["huge_formal_scale"] = huge_group is not None and huge_group["formal_powers"] == [10, -4, 0, 0]
    checks["member_refs_exact"] = all(
        group["member_product_group_id"] in product_groups
        and set(product_groups[group["member_product_group_id"]]["member_product_ids"]).issubset(set(result["products"]))
        for group in result["absolute_coefficient_bounds"]
    )

    class_only = _clone(packet)
    class_only["state_nodes"]["state:V4:W1:Z0:0"]["error_classes"] = ["fixture_v4:class-only"]
    class_only_result = build_fn(c, class_only, hb_admissibility=True, sc=sc, q_cells=q_cells)
    checks["state_error_class_without_flag_preserved"] = (
        class_only_result["products"]["v4p:positive"]["error_branch"] == "error"
        and "fixture_v4:class-only" in class_only_result["products"]["v4p:positive"]["state_error_classes"]
    )
    wrong_basis = _clone(packet)
    wrong_basis["source_power_basis"] = ["hb", "F0_anchor", "Lambda", "Pstar_squared"]
    bad_product_powers = _clone(packet)
    bad_product_powers["rows"]["h5"][0]["formal_powers"][0] += 1
    wrong_product_row = _clone(packet)
    wrong_product_row["rows"]["h5"][0]["row"] = "h5_Z1"
    checks["wrong_basis_rejected"] = _rejected(lambda: build_fn(c, wrong_basis, hb_admissibility=True, sc=sc))
    checks["wrong_product_powers_rejected"] = _rejected(lambda: build_fn(c, bad_product_powers, hb_admissibility=True, sc=sc))
    checks["wrong_product_row_rejected"] = _rejected(lambda: build_fn(c, wrong_product_row, hb_admissibility=True, sc=sc))

    missing_z1 = _clone(packet)
    del missing_z1["rows"]["h5_Z1"]
    missing_source = _clone(packet)
    del missing_source["source_nodes"]["v4node:positive"]
    unknown_family = _clone(packet)
    unknown_family["source_nodes"]["v4node:positive"]["source_family"] = "unknown"
    source_only_feedback = _clone(packet)
    source_only_feedback["rows"]["h5"][0]["source_only"] = True
    missing_state = _clone(packet)
    missing_state["rows"]["h5"][0]["state_node_id"] = None
    applied_scalar = _clone(packet)
    applied_scalar["rows"]["h5"][0]["scalar_applied_to_coefficient"] = True
    bad_exact = _clone(sc)
    del bad_exact["upper_exact_mpf_tuple"]
    bad_q = [
        _payload(c, 0, Fraction(1, 4)),
        _payload(c, Fraction(1, 4), Fraction(1, 2)),
        _payload(c, Fraction(1, 2), Fraction(3, 4)),
    ]
    checks["missing_Z1_rejected"] = _rejected(lambda: build_fn(c, missing_z1, hb_admissibility=True, sc=sc))
    checks["missing_source_rejected"] = _rejected(lambda: build_fn(c, missing_source, hb_admissibility=True, sc=sc))
    checks["unknown_family_rejected"] = _rejected(lambda: build_fn(c, unknown_family, hb_admissibility=True, sc=sc))
    checks["source_only_feedback_rejected"] = _rejected(lambda: build_fn(c, source_only_feedback, hb_admissibility=True, sc=sc))
    checks["missing_state_rejected"] = _rejected(lambda: build_fn(c, missing_state, hb_admissibility=True, sc=sc))
    checks["applied_scalar_rejected"] = _rejected(lambda: build_fn(c, applied_scalar, hb_admissibility=True, sc=sc))
    checks["missing_hb_certificate_rejected"] = _rejected(lambda: build_fn(c, packet, sc=sc))
    checks["false_hb_certificate_rejected"] = _rejected(lambda: build_fn(c, packet, hb_admissibility=False, sc=sc))
    checks["float_sc_rejected"] = _rejected(lambda: build_fn(c, packet, hb_admissibility=True, sc=1.0))
    checks["malformed_sc_rejected"] = _rejected(lambda: build_fn(c, packet, hb_admissibility=True, sc=bad_exact))
    checks["malformed_q_rejected"] = _rejected(lambda: build_fn(c, packet, hb_admissibility=True, sc=sc, q_cells=bad_q))

    scaled = build_fn(c, packet, hb_admissibility=True, sc=_payload(c, Fraction(1, 36), Fraction(1, 36)), q_cells=q_cells)
    base_S = _decode(c, result["cell_factors"]["q3"]["S_i"])
    scaled_S = _decode(c, scaled["cell_factors"]["q3"]["S_i"])
    checks["scale_changes_factor"] = scaled_S.b > base_S.b
    checks["scaled_keeps_four_cells"] = len(scaled["cell_factors"]) == 4

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "implementation_imported": False,
        "direct_K0_prefix_checked": True,
        "scalar_once_checked": True,
        "g5_axis_isolation_checked": True,
        "core_inlet_unknown_checked": True,
        "formal_hb_not_materialized_checked": True,
        "factorized_cell_groups_checked": True,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)


