"""Fresh symbolic P2 pressure-gap kernel draft.

This draft consumes only the saved P2 source-binding result. It does not
import, execute, copy, or replay any original provider, source, checker,
normalizer, reducer, quadrature, recurrence, primitive, sigma, root-shape, or
cocycle implementation. The proof objects are built from fresh small
polynomial operations over formal symbols. No source value is evaluated.
"""

from __future__ import annotations

import hashlib
from typing import Any, Callable, Mapping


INPUT_SCHEMA = "compact_O3_transition_P2_pressure_gap_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_P2_pressure_gap"
COMPONENT = "compact_O3_transition_P2_pressure_gap"

BINDING_RESULT_SHA256 = "007ffc9391381a6d62cfccb763dc6a4ce2b430cf51d6209013b576f2f36fc521"
BINDING_SOURCE_CHECKPOINT_SHA256 = "79a781ea0fa4a6509e907abb0f6b8c03f01f06ae1700bc9123e1f8064bf26a41"
BINDING_REPORT_SHA256 = "ac6768f0a842b8a0c3dba0fb8c24721bdce7c3c1435591c9bfdf6709db9dde3f"

ORDERS = (0, 1, 2, 3)
SIDES = ("native", "rebuilt_C3")
NEW_FLAGS = (
    "O3_transition_P2_signed_jet_gap_functions_defined",
    "O3_transition_P2_pressure_quotient_correlation_ordinary_Z3_proved",
)
GLOBAL_FALSE_FIELDS = (
    "global_five_moment_closure",
    "global_closure_proved",
    "O3_transition_P2_zero_claimed",
    "O3_transition_phase_equality_proved",
    "time_order_recursion_claimed",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "physical_layer_completed",
)
INHERITED_FALSE_FIELDS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
    "old_errors_removed",
    "old_error_terms_removed",
)
COUNT_KEYS = (
    "new_mathematical_proofs",
    "new_algebra_family_calls",
    "pressure_identity_orders",
    "gap_function_rows",
    "side_quotient_premise_rows",
    "numeric_calls",
    "source_provider_calls",
    "old_math_calls",
)

Poly = dict[tuple[str, ...], int]


class _Failure(Exception):
    pass


def _atom(name: str) -> Poly:
    return {(name,): 1}


def _zero() -> Poly:
    return {}


def _add(*parts: Poly) -> Poly:
    result: Poly = {}
    for part in parts:
        for monomial, coefficient in part.items():
            result[monomial] = result.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def _neg(part: Poly) -> Poly:
    return {monomial: -coefficient for monomial, coefficient in part.items()}


def _sub(left: Poly, right: Poly) -> Poly:
    return _add(left, _neg(right))


def _mul(left: Poly, right: Poly) -> Poly:
    if not left or not right:
        return {}
    result: Poly = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(sorted(left_monomial + right_monomial))
            result[monomial] = result.get(monomial, 0) + left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def _scale(part: Poly, coefficient: int) -> Poly:
    if coefficient == 0:
        return {}
    return {monomial: coefficient * value for monomial, value in part.items() if coefficient * value}


def _binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    value = 1
    for index in range(1, k + 1):
        value = value * (n - index + 1) // index
    return value


def _coefficient_family(left: list[Poly], right: list[Poly], order: int) -> Poly:
    return _add(
        *(
            _scale(_mul(left[index], right[order - index]), _binomial(order, index))
            for index in range(order + 1)
        )
    )


def _tree(part: Poly) -> dict[str, Any]:
    terms = []
    for monomial, coefficient in sorted(part.items(), key=lambda item: (len(item[0]), item[0])):
        terms.append(
            {
                "coefficient": coefficient,
                "factors": [{"kind": "formal_symbol", "name": name} for name in monomial],
            }
        )
    return {"kind": "coefficient_polynomial", "terms": terms, "is_zero": not terms}


def _tree_digest(part: Poly) -> str:
    return hashlib.sha256(repr(_tree(part)).encode("utf-8")).hexdigest()


def _symbol_tree(name: str) -> dict[str, str]:
    return {"kind": "formal_symbol", "name": name}


def _family_tree(left_name: str, right_name: str, order: int) -> dict[str, Any]:
    return {
        "kind": "ordinary_Z_binomial_product_family",
        "left_family": left_name,
        "right_family": right_name,
        "order": order,
        "terms": [
            {
                "left_order": index,
                "right_order": order - index,
                "binomial": _binomial(order, index),
            }
            for index in range(order + 1)
        ],
    }


def _recurrence_tree(order: int) -> dict[str, Any]:
    return {
        "kind": "ordinary_Z_signed_gap_source_function",
        "source_function": f"G_{order}",
        "order": order,
        "numerator": {
            "kind": "difference",
            "dW_source_function": f"dW_{order}",
            "dE_times_Pc_family": _family_tree("dE", "Pc", order),
            "previous_G_terms": [
                {
                    "binomial": _binomial(order, index),
                    "En_order": index,
                    "G_order": order - index,
                    "previous_G_source_function": f"G_{order - index}",
                }
                for index in range(1, order + 1)
            ],
        },
        "denominator": _symbol_tree("En_0"),
        "domain_premise": "En_0_nonzero",
    }


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise _Failure(f"missing mapping: {label}")
    return value


def _typed_view(row: Any, side: str, order: int) -> dict[str, Any]:
    view = _mapping(row, f"{side}.raw_W_derivative_views[{order}]")
    if view.get("order") != order or view.get("derivative_order") != order or view.get("dependency_order") != order + 1:
        raise _Failure(f"typed W source order changed: {side}:{order}")
    if view.get("p2_numerator_is_raw_W") is not False or view.get("derivative_creator_is_actual_source") is not True:
        raise _Failure(f"typed W source role changed: {side}:{order}")
    required = (
        "source_path",
        "source_sha",
        "source_function_ast_sha256",
        "source_assignment_ast_sha256",
    )
    if any(not isinstance(view.get(key), str) or not view.get(key) for key in required):
        raise _Failure(f"typed W source provenance missing: {side}:{order}")
    return {
        "side": side,
        "order": order,
        "dependency_order": order + 1,
        "source_path": view["source_path"],
        "source_sha": view["source_sha"],
        "source_function_ast_sha256": view["source_function_ast_sha256"],
        "source_assignment_ast_sha256": view["source_assignment_ast_sha256"],
        "p2_numerator_is_raw_W": False,
        "derivative_creator_is_actual_source": True,
        "actual_value_present": False,
    }


def _compact_literal_handles(side: str, reference: Mapping[str, Any]) -> dict[str, Any]:
    literal = _mapping(reference.get("literal_handles"), f"{side}.literal_handles")
    selectors = literal.get("p2_root_selectors")
    if not isinstance(selectors, list) or len(selectors) != len(ORDERS):
        raise _Failure(f"P2 literal selector count changed: {side}")
    compact_selectors = []
    for order, selector in enumerate(selectors):
        row = _mapping(selector, f"{side}.p2_root_selectors[{order}]")
        if row.get("order") != order or row.get("dependency_order") != order + 1:
            raise _Failure(f"P2 literal selector order changed: {side}:{order}")
        if row.get("root_numerator_is_raw_W") is not False:
            raise _Failure(f"P2 quotient selector was treated as W: {side}:{order}")
        compact_selectors.append(
            {
                key: row[key]
                for key in (
                    "order",
                    "dependency_order",
                    "root_node",
                    "recipe_node",
                    "root_numerator_is_raw_W",
                    "p2_numerator_role",
                    "quantity_path",
                    "Taylor_coefficient_factorial",
                    "factorial_applied_once",
                )
                if key in row
            }
        )
    return {
        "side": side,
        "archive_path": literal.get("archive_path"),
        "archive_sha256": literal.get("archive_sha256"),
        "chart": literal.get("chart"),
        "coordinate_node": literal.get("coordinate_node"),
        "Z_variable_node": literal.get("Z_variable_node"),
        "root_ids": list(literal.get("root_ids", [])),
        "positive_E_denominator_node": literal.get("positive_E_denominator_node"),
        "factorial_map": literal.get("factorial_map"),
        "p2_root_selectors": compact_selectors,
    }


def _compact_formula_route(side: str, reference: Mapping[str, Any]) -> dict[str, Any]:
    formula = _mapping(reference.get("formula"), f"{side}.formula")
    if side == "native":
        required = ("source_path", "source_sha")
        if any(not isinstance(formula.get(key), str) or not formula.get(key) for key in required):
            raise _Failure(f"source W formula route missing: {side}")
        return {
            "W_source_path": formula["source_path"],
            "W_source_sha": formula["source_sha"],
            "W_function": formula.get("function"),
            "W_function_ast_sha256": formula.get("function_ast_sha256"),
            "E_route_kind": "positive_E_denominator_source_handle",
            "positive_E_denominator_node": reference.get("literal_handles", {}).get("positive_E_denominator_node"),
            "actual_value_present": False,
        }
    raw_W = _mapping(formula.get("raw_W_source"), f"{side}.raw_W_source")
    if any(not isinstance(raw_W.get(key), str) or not raw_W.get(key) for key in ("source_path", "source_sha")):
        raise _Failure(f"source W formula route missing: {side}")
    quotient_route = _mapping(formula.get("source_quotients_route"), f"{side}.source_quotients_route")
    E_ast = quotient_route.get("E_quotient_call_ast_sha256")
    if not isinstance(E_ast, str) or not E_ast:
        raise _Failure("rebuilt_C3 E quotient source route missing")
    if any(not isinstance(quotient_route.get(key), str) or not quotient_route.get(key) for key in ("source_path", "source_sha")):
        raise _Failure("rebuilt_C3 E quotient source path missing")
    return {
        "W_source_path": raw_W["source_path"],
        "W_source_sha": raw_W["source_sha"],
        "W_function": raw_W.get("function"),
        "W_function_ast_sha256": raw_W.get("function_ast_sha256"),
        "E_route_kind": "source_quotients_E_quotient_call",
        "E_source_path": quotient_route["source_path"],
        "E_source_sha": quotient_route["source_sha"],
        "E_quotient_call_ast_sha256": E_ast,
        "source_quotients_function_ast_sha256": quotient_route.get("function_ast_sha256"),
        "actual_value_present": False,
    }


def _validated_binding(source: Mapping[str, Any]) -> dict[str, Any]:
    expected_pins = {
        "binding_result_raw_sha256": BINDING_RESULT_SHA256,
        "binding_source_checkpoint_sha256": BINDING_SOURCE_CHECKPOINT_SHA256,
        "binding_report_raw_sha256": BINDING_REPORT_SHA256,
    }
    for key, expected in expected_pins.items():
        if source.get(key) != expected:
            raise _Failure(f"binding raw SHA pin changed: {key}")
    binding = _mapping(source.get("binding_result"), "binding_result")
    flags = binding.get("new_flags")
    if not isinstance(flags, Mapping):
        flags = binding.get("flags")
    required_flags = (
        "O3_transition_P2_source_ports_bound",
        "O3_transition_P2_pressure_restoration_route_bound",
        "O3_transition_P2_ordinary_Z_source_projection_bound",
    )
    if not isinstance(flags, Mapping) or any(flags.get(key) is not True for key in required_flags):
        raise _Failure("saved binding flags are not all true")
    references = _mapping(binding.get("typed_source_references"), "binding_result.typed_source_references")
    views = _mapping(binding.get("whole_function_pressure_numerator_views"), "binding_result.whole_function_pressure_numerator_views")
    if set(references) != set(SIDES) or set(views) != set(SIDES):
        raise _Failure("saved binding sides are not exactly native and rebuilt_C3")
    typed_W: dict[str, list[dict[str, Any]]] = {}
    literals: dict[str, dict[str, Any]] = {}
    E_routes: dict[str, dict[str, Any]] = {}
    for side in SIDES:
        reference = _mapping(references[side], f"{side}.typed_reference")
        side_views = _mapping(views[side], f"{side}.whole_function_views")
        rows = side_views.get("raw_W_derivative_views")
        if not isinstance(rows, list) or len(rows) != len(ORDERS):
            raise _Failure(f"saved binding W view count changed: {side}")
        typed_W[side] = [_typed_view(row, side, order) for order, row in enumerate(rows)]
        literals[side] = _compact_literal_handles(side, reference)
        E_routes[side] = _compact_formula_route(side, reference)
    return {
        "binding_result_raw_sha256": BINDING_RESULT_SHA256,
        "binding_source_checkpoint_sha256": BINDING_SOURCE_CHECKPOINT_SHA256,
        "binding_report_raw_sha256": BINDING_REPORT_SHA256,
        "typed_W_views": typed_W,
        "P2_literal_handles": literals,
        "source_E_routes": E_routes,
        "binding_flags": {key: True for key in required_flags},
    }


def _source_E_views(
    typed_W: Mapping[str, Any],
    E_routes: Mapping[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    for side in SIDES:
        result[side] = []
        for order, W_view in enumerate(typed_W[side]):
            result[side].append(
                {
                    "side": side,
                    "order": order,
                    "dependency_order": order + 1,
                    "symbol": f"E{'n' if side == 'native' else 'c'}_{order}",
                    "source_expression_kind": "typed_source_E_whole_function_derivative_expression",
                    "source_expression_label": "sourceEn" if side == "native" else "sourceEc",
                    "source_E_route": E_routes[side],
                    "W_view_route": {
                        "source_path": W_view["source_path"],
                        "source_sha": W_view["source_sha"],
                        "source_function_ast_sha256": W_view["source_function_ast_sha256"],
                        "source_assignment_ast_sha256": W_view["source_assignment_ast_sha256"],
                    },
                    "source_zero_typed": order == 0,
                    "actual_value_present": False,
                }
            )
    return result


def _side_symbols(side: str, quantity: str) -> list[Poly]:
    prefix = {"native": "n", "rebuilt_C3": "c"}[side]
    return [_atom(f"{quantity}{prefix}_{order}") for order in ORDERS]


def _side_quotient_premises(
    typed_W: Mapping[str, Any],
    literals: Mapping[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    for side in SIDES:
        E = _side_symbols(side, "E")
        P = _side_symbols(side, "P")
        W = _side_symbols(side, "W")
        for order in ORDERS:
            rhs = _coefficient_family(E, P, order)
            rows.append(
                {
                    "side": side,
                    "order": order,
                    "dependency_order": order + 1,
                    "source_W_symbol": _tree(W[order]),
                    "quotient_premise_rhs": _tree(rhs),
                    "quotient_premise_family": _family_tree(
                        f"E{'n' if side == 'native' else 'c'}",
                        f"P{'n' if side == 'native' else 'c'}",
                        order,
                    ),
                    "source_W_view": typed_W[side][order],
                    "P2_literal_selector": literals[side]["p2_root_selectors"][order],
                    "premise_used_for_substitution_only": True,
                    "W_source_producer_relation_is_not_extra_independent_physics": True,
                    "W_is_not_defined_as_E_times_P2": True,
                }
            )
    return rows


def _gap_rows(typed_W: Mapping[str, Any]):
    Wn = _side_symbols("native", "W")
    Wc = _side_symbols("rebuilt_C3", "W")
    En = _side_symbols("native", "E")
    Ec = _side_symbols("rebuilt_C3", "E")
    Pn = _side_symbols("native", "P")
    Pc = _side_symbols("rebuilt_C3", "P")
    dW = [_sub(Wn[order], Wc[order]) for order in ORDERS]
    dE = [_sub(En[order], Ec[order]) for order in ORDERS]
    for order in ORDERS:
        independent_G = _sub(Pn[order], Pc[order])
        G_witness = [_sub(Pn[index], Pc[index]) for index in ORDERS]
        lhs = _coefficient_family(En, G_witness, order)
        rhs_before = _sub(dW[order], _coefficient_family(dE, Pc, order))
        Wn_quotient = _coefficient_family(En, Pn, order)
        Wc_quotient = _coefficient_family(Ec, Pc, order)
        dW_after = _sub(Wn_quotient, Wc_quotient)
        rhs_after = _sub(dW_after, _coefficient_family(dE, Pc, order))
        residual = _sub(lhs, rhs_after)
        if residual:
            raise _Failure(f"fresh pressure quotient residual did not cancel at order {order}")
        recurrence_numerator = _sub(
            dW[order],
            _add(
                _coefficient_family(dE, Pc, order),
                *(
                    _scale(_mul(En[index], _atom(f"G_{order - index}")), _binomial(order, index))
                    for index in range(1, order + 1)
                ),
            ),
        )
        recurrence_numerator_after_substitution = _sub(
            dW_after,
            _add(
                _coefficient_family(dE, Pc, order),
                *(
                    _scale(_mul(En[index], G_witness[order - index]), _binomial(order, index))
                    for index in range(1, order + 1)
                ),
            ),
        )
        induction_target = _mul(En[0], independent_G)
        induction_residual = _sub(recurrence_numerator_after_substitution, induction_target)
        if induction_residual:
            raise _Failure(f"fresh pressure quotient induction residual did not cancel at order {order}")
        previous = [
            {
                "order": order - index,
                "source_function": f"G_{order - index}",
                "binomial": _binomial(order, index),
            }
            for index in range(1, order + 1)
        ]
        yield {
            "order": order,
            "dependency_order": order + 1,
            "source_function": f"G_{order}",
            "relation": "En*G=dW-dE*Pc",
            "domain_premises": [
                {"kind": "nonzero", "symbol": "En_0", "side": "native", "actual_value_present": False},
                {"kind": "nonzero", "symbol": "Ec_0", "side": "rebuilt_C3", "actual_value_present": False},
            ],
            "same_Z_chart_pairing": {
                "formal_pairing": True,
                "physical_coordinate_equality_asserted": False,
                "phase_equality_asserted": False,
            },
            "independent_gap_witness": _tree(independent_G),
            "lhs_derivative_construct": {
                "family": _family_tree("En", "G", order),
                "coefficient_polynomial": _tree(lhs),
            },
            "rhs_derivative_construct_before_substitution": {
                "dW_coefficient": _tree(dW[order]),
                "dE_times_Pc_family": _family_tree("dE", "Pc", order),
                "coefficient_polynomial": _tree(rhs_before),
            },
            "source_quotient_premise_substitution": {
                "native_W_coefficient": _tree(Wn_quotient),
                "rebuilt_C3_W_coefficient": _tree(Wc_quotient),
                "dW_coefficient_after_substitution": _tree(dW_after),
                "premise_is_not_new_W_definition": True,
                "source_producer_relation_is_conditional": True,
            },
            "rhs_derivative_construct_after_substitution": {
                "coefficient_polynomial": _tree(rhs_after),
            },
            "residual_polynomial": _tree(residual),
            "residual_digest": _tree_digest(residual),
            "residual_is_zero": True,
            "source_recurrence_gap_function": _recurrence_tree(order),
            "source_recurrence_numerator_expanded_with_previous_G_atoms": _tree(recurrence_numerator),
            "recurrence_induction_witness": {
                "premise_substituted_previous_G": [
                    {
                        "order": order - index,
                        "source_function": f"G_{order - index}",
                        "independent_witness": _tree(G_witness[order - index]),
                        "binomial": _binomial(order, index),
                    }
                    for index in range(1, order + 1)
                ],
                "left_recurrence_numerator_after_quotient_and_previous_G_substitution": _tree(recurrence_numerator_after_substitution),
                "right_En0_times_independent_gap_witness": _tree(induction_target),
                "residual_polynomial": _tree(induction_residual),
                "residual_digest": _tree_digest(induction_residual),
                "residual_is_zero": True,
                "induction_is_separate_from_recurrence_self_subtraction": True,
            },
            "previous_G_source_functions": previous,
            "independent_witness_not_recurrence_self_subtraction": True,
            "typed_W_views": {
                "native": typed_W["native"][order],
                "rebuilt_C3": typed_W["rebuilt_C3"][order],
            },
        }


def _counts(
    accepted: bool,
    completed_orders: int = 0,
    premise_rows: int = 0,
    algebra_entered: bool = False,
) -> dict[str, int]:
    result = {key: 0 for key in COUNT_KEYS}
    if accepted:
        result.update(
            new_mathematical_proofs=1,
            new_algebra_family_calls=1,
            pressure_identity_orders=4,
            gap_function_rows=4,
            side_quotient_premise_rows=8,
        )
    elif algebra_entered or completed_orders or premise_rows:
        result.update(
            new_algebra_family_calls=1,
            pressure_identity_orders=completed_orders,
            gap_function_rows=completed_orders,
            side_quotient_premise_rows=premise_rows,
        )
    return result


def _audit_zero() -> dict[str, int]:
    return {
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "provider_calls": 0,
        "checker_calls": 0,
        "normalizer_calls": 0,
        "reducer_calls": 0,
        "quadrature_calls": 0,
        "recurrence_calls": 0,
        "field_evaluations": 0,
        "integral_evaluations": 0,
    }


def _failure(
    events: list[dict[str, Any]],
    error: str,
    source: Mapping[str, Any] | None,
    completed_orders: int = 0,
    premise_rows: int = 0,
    algebra_entered: bool = False,
    completed_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    flags = {key: False for key in NEW_FLAGS}
    counts = _counts(False, completed_orders, premise_rows, algebra_entered)
    result: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "failed",
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "runtime_model_metadata": None,
        "events": list(events),
        "journal_events_prefix": list(events),
        "failure_prefix_events": list(events),
        "completed_proof_rows": list(completed_rows or []),
        "completed_recurrence_induction_rows": [
            row["recurrence_induction_witness"]
            for row in (completed_rows or [])
            if isinstance(row, Mapping) and "recurrence_induction_witness" in row
        ],
        "blockers": [error],
        "flags": dict(flags),
        "new_flags": dict(flags),
        **flags,
        "actual_counts": counts,
        "audit_counts": dict(counts),
        "new_mathematical_proofs": counts["new_mathematical_proofs"],
        "new_algebra_family_calls": counts["new_algebra_family_calls"],
        "pressure_identity_orders": counts["pressure_identity_orders"],
        "gap_function_rows": counts["gap_function_rows"],
        "side_quotient_premise_rows": counts["side_quotient_premise_rows"],
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "audit": _audit_zero(),
        "source_packet_defined": False,
        "source_packet_definition_accepted": False,
        "source_packet_accepted_owned": False,
        "source_packet_accepted_inherited": False,
        "full_goal_complete": False,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "inherited_global_flags_all_false": True,
        "inherited_false_flags": {key: False for key in (*GLOBAL_FALSE_FIELDS, *INHERITED_FALSE_FIELDS)},
        **{key: False for key in GLOBAL_FALSE_FIELDS},
        **{key: False for key in INHERITED_FALSE_FIELDS},
    }
    if isinstance(source, Mapping):
        for key in (
            "binding_result_raw_sha256",
            "binding_source_checkpoint_sha256",
            "binding_report_raw_sha256",
        ):
            result[key] = source.get(key)
    return result


def derive(
    source: Mapping[str, Any],
    diagnostic: Callable[[dict[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Build four formal same-Z pressure-gap jet rows from saved source ports."""

    events: list[dict[str, Any]] = []
    completed_orders = 0
    premise_rows = 0
    algebra_entered = False
    completed_rows: list[dict[str, Any]] = []

    def emit(event: dict[str, Any]) -> None:
        events.append(event)
        if diagnostic is not None:
            if not callable(diagnostic):
                raise _Failure("diagnostic callback must be callable")
            diagnostic(event)

    try:
        if not isinstance(source, Mapping) or source.get("schema") != INPUT_SCHEMA:
            raise _Failure("fresh pressure-gap API schema changed")
        if diagnostic is not None and not callable(diagnostic):
            raise _Failure("diagnostic callback must be callable")
        binding = _validated_binding(source)
        emit(
            {
                "event": "saved_P2_binding_validated",
                "typed_W_view_count": 8,
                "binding_flags": binding["binding_flags"],
                "binding_result_raw_sha256": binding["binding_result_raw_sha256"],
            }
        )
        typed_W = binding["typed_W_views"]
        source_E = _source_E_views(typed_W, binding["source_E_routes"])
        algebra_entered = True
        emit(
            {
                "event": "fresh_pressure_quotient_algebra_entered",
                "ordinary_Z_orders": list(ORDERS),
                "side_quotient_premise_rows_planned": 8,
                "side_quotient_premise_rows_actual": None,
            }
        )
        premises = _side_quotient_premises(typed_W, binding["P2_literal_handles"])
        premise_rows = len(premises)
        gap_rows = []
        for row in _gap_rows(typed_W):
            gap_rows.append(row)
            completed_rows.append(row)
            completed_orders += 1
            emit(
                {
                    "event": "pressure_gap_proof_row",
                    "order": row["order"],
                    "dependency_order": row["dependency_order"],
                    "source_function": row["source_function"],
                    "residual_is_zero": row["residual_is_zero"],
                    "induction_residual_is_zero": row["recurrence_induction_witness"]["residual_is_zero"],
                    "residual_digest": row["residual_digest"],
                    "previous_G_source_functions": row["previous_G_source_functions"],
                    "proof_row": row,
                }
            )
        flags = {key: True for key in NEW_FLAGS}
        counts = _counts(True)
        audit = _audit_zero()
        return {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "status": "passed",
            "passed": True,
            "source_only": True,
            "metadata_only": False,
            "runtime_model_metadata": None,
            "events": events,
            "journal_events": list(events),
            "flags": dict(flags),
            "new_flags": dict(flags),
            **flags,
            "actual_counts": counts,
            "audit_counts": dict(counts),
            "new_mathematical_proofs": counts["new_mathematical_proofs"],
            "new_algebra_family_calls": counts["new_algebra_family_calls"],
            "pressure_identity_orders": counts["pressure_identity_orders"],
            "gap_function_rows": counts["gap_function_rows"],
            "side_quotient_premise_rows": counts["side_quotient_premise_rows"],
            "numeric_calls": 0,
            "source_provider_calls": 0,
            "old_math_calls": 0,
            "audit": audit,
            "binding_pins": {
                "binding_result_raw_sha256": binding["binding_result_raw_sha256"],
                "binding_source_checkpoint_sha256": binding["binding_source_checkpoint_sha256"],
                "binding_report_raw_sha256": binding["binding_report_raw_sha256"],
            },
            "typed_source_W_views": typed_W,
            "actual_W_source_producers": {
                side: {
                    "construction_kind": "typed_complete_fullpressure_source_producer",
                    "source_assignment_views": typed_W[side],
                    "not_recovered_from_E_times_P2": True,
                }
                for side in SIDES
            },
            "typed_source_E_views": source_E,
            "formal_side_quotients": {
                "native": {
                    "P_symbol_prefix": "Pn",
                    "W_symbol_prefix": "Wn",
                    "E_symbol_prefix": "En",
                    "relation_kind": "formal_source_quotient_Pn_equals_Wn_over_En",
                    "actual_value_present": False,
                    "source_W_assignment_independent": True,
                },
                "rebuilt_C3": {
                    "P_symbol_prefix": "Pc",
                    "W_symbol_prefix": "Wc",
                    "E_symbol_prefix": "Ec",
                    "relation_kind": "formal_source_quotient_Pc_equals_Wc_over_Ec",
                    "actual_value_present": False,
                    "source_W_assignment_independent": True,
                },
            },
            "P2_literal_handles": binding["P2_literal_handles"],
            "side_quotient_premises": premises,
            "pressure_identity_rows": gap_rows,
            "gap_function_rows_detail": gap_rows,
            "recurrence_induction_rows": [row["recurrence_induction_witness"] for row in gap_rows],
            "ordinary_Z_orders": list(ORDERS),
            "domain_premises": [
                {"kind": "nonzero", "symbol": "En_0", "side": "native", "actual_value_present": False},
                {"kind": "nonzero", "symbol": "Ec_0", "side": "rebuilt_C3", "actual_value_present": False},
            ],
            "same_Z_chart_pairing": {
                "formal_pairing": True,
                "physical_radius_equality_asserted": False,
                "phase_equality_asserted": False,
                "geometry_history_pressure_independence_preserved": True,
            },
            "pressure_identity": "En*G=dW-dE*Pc",
            "source_quotient_substitution_is_premise_not_definition": True,
            "source_quotient_relation_is_conditional_on_side_source": True,
            "independent_gap_witness_is_not_recurrence_selfzero": True,
            "new_algebra_family": "fresh_symbolic_pressure_quotient_correlation",
            "source_packet_defined": True,
            "source_packet_definition_accepted": True,
            "source_packet_accepted_owned": True,
            "source_packet_accepted_inherited": False,
            "full_goal_complete": False,
            "numeric_values_computed": False,
            "numeric_field_values_present": False,
            "actual_numeric_jets_available": False,
            "inherited_global_flags_all_false": True,
            "inherited_false_flags": {key: False for key in (*GLOBAL_FALSE_FIELDS, *INHERITED_FALSE_FIELDS)},
            **{key: False for key in GLOBAL_FALSE_FIELDS},
            **{key: False for key in INHERITED_FALSE_FIELDS},
            "unresolved": [
                "No P2 value, pressure closure, phase equality, radius equality, or signed-gap zero is asserted.",
                "The quotient identities are formal same-Z premises; no original source provider is evaluated.",
            ],
        }
    except _Failure as exc:
        return _failure(
            events,
            str(exc),
            source if isinstance(source, Mapping) else None,
            completed_orders,
            premise_rows,
            algebra_entered,
            completed_rows,
        )
    except Exception as exc:
        return _failure(
            events,
            f"unexpected fresh pressure-gap failure: {type(exc).__name__}: {exc}",
            source if isinstance(source, Mapping) else None,
            completed_orders,
            premise_rows,
            algebra_entered,
            completed_rows,
        )



