"""Fresh, source-only first radial variation compiler for the O3 chart.

This module deliberately builds a symbolic value/tangent IR.  It never imports
the project source, evaluates a source expression, materializes ``N``, or
turns an enclosure into a value.  ``derive`` consumes the literal capsule
whose graph records are supplied by the caller and returns a serializable
report suitable for a later, separately controlled invocation.

The accepted capsule shape is intentionally small and explicit::

    {
        "source_graph_records": {"node": {"operation": ...}},
        "roots": {"E_N": node, "V_N": node, ...},
        "frame": {"N": 18, "native_x": 241536, ...},
        "held_source_input_nodes": [...],
        "source_history_tangents": {"node": {"source_expression": ...}},
        "source_premises": [...],
    }

The source graph is a bounded literal DAG.  Prefix nodes below the frame's
``fixed_prefix_upper_exclusive`` are retained as opaque source references.
Their local source-domain tangent is a symbolic data reference, never a
numeric zero and never a claim of global incoming-x independence.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_radial_variation_kernel_result_v1"
FAMILY = "native_O3_radial_first_variation"
COMPONENT = "native_O3_radial_variation"

# These are deliberately literal audit contracts.  The durable runner reads
# them as metadata without importing or executing this compiler.  The first
# three family counters are reached by the family-entered/completed callback
# rows; the remaining zero counters retain the no-numeric/no-ancestor policy.
SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "actual_negative_api_refusals": 0,
    "cached_algebra_replay_calls": 0,
    "failed_family_calls": 0,
    "new_algebra_family_calls": 1,
    "new_mathematical_proofs": 1,
    "new_source_identity_family_calls": 1,
    "numeric_evaluation_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "source_family_calls": 1,
    "source_value_identity_calls": 0,
}

SUCCESS_FLAGS = {
    "actual_E_background_exp_J_sigma_tangent": True,
    "actual_A_tangent": True,
    "actual_B_over_Pstar_tangent": True,
    "actual_E_N_tangent": True,
    "actual_V_N_tangent": True,
    "actual_phase_fractional_part_tangent": True,
    "conditional_C2_regular_open_cell": True,
    "conditional_native_O3_first_variation_compiled": True,
    "conditional_positive_inverse_derivative": True,
    "actual_incoming_x_independence_proved": False,
    "global_phase_seam_regularity_proved": False,
    "global_axis_regularity_proved": False,
    "global_inverse_conditioning_proved": False,
    "physical_B_confused_with_B_over_Pstar": False,
    "N_materialized": False,
    "numeric_evaluation": False,
    "oracle_source_value_identity": False,
    "native_field_derivative_identity_proved": False,
    "cap_deletion": False,
    "temporal_recursion": False,
    "physical_layer_completed": False,
    "fullphysical_construction": False,
    "full_goal_complete": False,
}

CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "source_history_join_summary": 1,
    "phase_variation_summary": 1,
    "branch_integral_sigma_summary": 1,
    "inverse_ift_summary": 1,
    "finite_N_definition": 1,
    "family_completed": 1,
}

COUNT_EVENT_MAP = {
    "family_entered": [
        "new_source_identity_family_calls",
        "new_algebra_family_calls",
        "source_family_calls",
    ],
    "source_history_join_summary": [],
    "phase_variation_summary": [],
    "branch_integral_sigma_summary": [],
    "inverse_ift_summary": [],
    "finite_N_definition": [],
    "family_completed": ["new_mathematical_proofs"],
}

SUPPORTED_SOURCE_OPERATIONS = frozenset(
    {
        "sum",
        "product",
        "positive_quotient",
        "negative",
        "analytic_unary",
        "exact_original_source_flat_sigma",
        "bound_variable",
        "original_lazy_flat_branch",
        "substitute_original_inverse_angle",
        "original_monotone_phase_inverse",
        "definite_integral",
        "function_substitution",
        "exact_real_comparison",
        "exact_rational",
        "exact_positive_integer_expression",
        "current_original_source_parameter",
        "current_original_leading_function_recipe",
    }
)

ANALYTIC_NAMES = frozenset({"exp", "log", "positive_sqrt", "cos", "sin", "fractional_part"})


class KernelInputError(ValueError):
    """Raised internally when a literal graph violates the narrow contract."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _node_number(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("-"):
            body = text[1:]
        else:
            body = text
        if body.isdigit() and body:
            try:
                return int(text)
            except ValueError:
                return None
    return None


def _scalar(value: Any) -> Any:
    """Return a JSON-safe scalar for compact receipts without coercing math."""

    if value is None or isinstance(value, (str, int, bool, float)):
        return value
    return str(value)


class _Context:
    __slots__ = ("mode", "held_angle", "dummy_variables")

    def __init__(
        self,
        mode: str,
        held_angle: str | None = None,
        dummy_variables: tuple[str, ...] = (),
    ) -> None:
        if mode not in {"y", "angle"}:
            raise KernelInputError("unknown tangent mode " + repr(mode))
        self.mode = mode
        self.held_angle = held_angle
        self.dummy_variables = tuple(dummy_variables)

    def with_held_angle(self, name: str | None) -> "_Context":
        return _Context(self.mode, name, self.dummy_variables)

    def with_dummy(self, name: str) -> "_Context":
        if name in self.dummy_variables:
            return self
        return _Context(self.mode, self.held_angle, self.dummy_variables + (name,))

    def key(self) -> tuple[str, str | None, tuple[str, ...]]:
        return self.mode, self.held_angle, self.dummy_variables


class _Compiler:
    """Build the value and first-variation expression graphs."""

    def __init__(
        self,
        capsule: Mapping[str, Any],
        saved_flux_result: Mapping[str, Any] | None,
        callback: Callable[[str, Mapping[str, Any]], Any] | None,
    ) -> None:
        self.capsule = _as_mapping(capsule, "capsule")
        if saved_flux_result is not None and not _is_mapping(saved_flux_result):
            raise KernelInputError("saved_flux_result must be a mapping or None")
        self.saved = saved_flux_result if saved_flux_result is not None else {}
        if callback is not None and not callable(callback):
            raise KernelInputError("callback must be callable or None")
        self.callback = callback

        records = self.capsule.get("source_graph_records")
        if not _is_mapping(records) or not records:
            raise KernelInputError("source_graph_records must be a nonempty mapping")
        self.records: Mapping[Any, Any] = records
        self.record_keys = {str(key) for key in records}

        self.frame = _as_mapping(self.capsule.get("frame"), "frame")
        self.roots = _as_mapping(self.capsule.get("roots"), "roots")
        self.history = _as_mapping(
            self.capsule.get("source_history_tangents", {}), "source_history_tangents"
        )
        self.held_inputs = {
            self._id_string(item)
            for item in self.capsule.get("held_source_input_nodes", ())
        }
        self.prefix_upper = self._required_frame_int(
            "fixed_prefix_upper_exclusive", 241536
        )
        self.x_node = self._required_frame_int("native_x", 241536)
        self.angle_node = self._required_frame_int("angle", 243033)
        self.phase_node = self._required_frame_int("actual_phase", 242551)
        self.N_node = self._required_frame_int("N", 18)
        self.R_node = self._required_frame_int("native_R", 241787)
        self.angle_name = self._bound_name(self.angle_node)

        self.value_ir: list[dict[str, Any]] = []
        self.tangent_ir: list[dict[str, Any]] = []
        self.equations: list[dict[str, Any]] = []
        self.events: list[dict[str, Any]] = []
        self.value_cache: dict[tuple[Any, ...], dict[str, Any]] = {}
        self.tangent_cache: dict[tuple[Any, ...], dict[str, Any]] = {}
        self.source_history_refs: dict[str, dict[str, Any]] = {}
        self.touched_nodes: set[str] = set()
        self.opaque_inputs: set[str] = set()
        self._value_active: set[tuple[Any, ...]] = set()
        self._tangent_active: set[tuple[Any, ...]] = set()
        self._value_counter = 0
        self._tangent_counter = 0
        self.counts: dict[str, int] = {
            "source_graph_nodes_seen": 0,
            "value_ir_rows": 0,
            "tangent_ir_rows": 0,
            "equation_receipts": 0,
            "source_history_joins": 0,
            "prefix_source_references": 0,
            "inverse_ift_rows": 0,
            "integral_leibniz_rows": 0,
            "lazy_branch_rows": 0,
            "sigma_variation_rows": 0,
            "phase_fractional_rows": 0,
            "unknown_operation_refusals": 0,
            "callback_events": 0,
        }
        self.flags: dict[str, Any] = {
            "conditional_native_O3_first_variation_compiled": False,
            "conditional_C2_regular_open_cell": False,
            "conditional_positive_inverse_derivative": False,
            "actual_E_background_exp_J_sigma_tangent": False,
            "actual_A_tangent": False,
            "actual_B_over_Pstar_tangent": False,
            "actual_E_N_tangent": False,
            "actual_V_N_tangent": False,
            "actual_phase_fractional_part_tangent": False,
            "native_y_equals_D_x": bool(self.frame.get("total_y_equals_D_x") is True),
            "N_materialized": False,
            "numeric_evaluation": False,
            "actual_incoming_x_independence_proved": False,
            "global_phase_seam_regularity_proved": False,
            "global_axis_regularity_proved": False,
            "global_inverse_conditioning_proved": False,
            "physical_B_confused_with_B_over_Pstar": False,
            "cap_deletion": False,
            "temporal_recursion": False,
            "fullphysical_construction": False,
            "full_goal_complete": False,
        }
        self.premises = list(self.capsule.get("source_premises", ()))
        self._validate_frame()
        self._emit(
            "family_entered",
            {
                "family": FAMILY,
                "schema": self.capsule.get("schema"),
                "source_only": True,
                "numeric_evaluation": False,
                "N_materialized": False,
                "regularity_scope": "conditional C2 regular open O3 phase cell",
                "x_scope": "0<x<1",
                "phase_scope": "actual fractional_part source phase off seams",
                "inverse_scope": "positive Phi_angle source premise",
                "ordinary_Z": self.frame.get("ordinary_Z"),
            },
        )

    @staticmethod
    def _id_string(value: Any) -> str:
        node = _node_number(value)
        return str(node) if node is not None else str(value)

    def _required_frame_int(self, key: str, default: int) -> int:
        raw = self.frame.get(key, default)
        value = _node_number(raw)
        if value is None:
            raise KernelInputError("frame." + key + " must be an opaque node reference")
        return value

    def _validate_frame(self) -> None:
        if self.prefix_upper <= 0:
            raise KernelInputError("fixed prefix boundary must be positive")
        for key, expected in (
            ("N", 18),
            ("native_x", 241536),
            ("angle", 243033),
            ("actual_phase", 242551),
        ):
            value = _node_number(self.frame.get(key, expected))
            if value != expected:
                raise KernelInputError("frame." + key + " changed from required literal")
        held_numbers = {_node_number(item) for item in self.capsule.get("held_source_input_nodes", ())}
        if self.N_node not in held_numbers:
            raise KernelInputError("opaque N reference is not retained as a held source input")
        if self.frame.get("physical_theta_is_distinct") is not True:
            raise KernelInputError("physical theta must remain distinct from phase phi")
        if self.frame.get("total_y_equals_D_x") is not True:
            raise KernelInputError("native y=D_x premise is required")

    def _raw_row(self, ref: Any) -> Mapping[str, Any] | None:
        if _is_mapping(ref):
            if "operation" in ref or "op" in ref:
                return ref
            for key in ("node", "node_id", "id", "source_node", "ref"):
                if key in ref:
                    return self._raw_row(ref[key])
            return None
        key = str(ref)
        if key not in self.record_keys:
            return None
        value = self.records.get(ref)
        if value is None and key != ref:
            value = self.records.get(key)
        return value if _is_mapping(value) else None

    def _ref_id(self, ref: Any) -> int | None:
        if _is_mapping(ref):
            for key in ("node", "node_id", "id", "source_node", "ref"):
                if key in ref:
                    result = _node_number(ref[key])
                    if result is not None:
                        return result
            return None
        return _node_number(ref)

    def _ref_key(self, ref: Any) -> str:
        node = self._ref_id(ref)
        if node is not None:
            return str(node)
        if _is_mapping(ref):
            return "inline:" + str(sorted((str(k), str(v)) for k, v in ref.items()))
        return str(ref)

    def _bound_name(self, ref: Any) -> str | None:
        row = self._raw_row(ref)
        if row is not None and row.get("operation") == "bound_variable":
            name = row.get("name")
            return str(name) if name is not None else None
        if isinstance(ref, str):
            return ref
        return None

    def _op(self, row: Mapping[str, Any], source_node: Any) -> str:
        operation = row.get("operation", row.get("op"))
        if not isinstance(operation, str):
            self.counts["unknown_operation_refusals"] += 1
            raise KernelInputError("source node " + str(source_node) + " has no operation")
        if operation not in SUPPORTED_SOURCE_OPERATIONS:
            self.counts["unknown_operation_refusals"] += 1
            raise KernelInputError(
                "unknown source operation " + repr(operation) + " at node " + str(source_node)
            )
        return operation

    def _provenance(self, row: Mapping[str, Any]) -> dict[str, Any]:
        keep = (
            "name",
            "definition",
            "source_binding",
            "source_positive_certificate",
            "derivative_positive_certificate",
            "measure",
            "variable",
            "ordinary_Z_order",
            "Z_order",
            "native_chart",
            "exact_function_integral",
            "exact_original_sigma_integral",
            "angle_endpoints_independent_of_Z",
            "active_body_not_evaluated_on_flat_branch",
            "same_original_flat_extension_through_order5",
            "exact_function_not_its_interval_cap",
        )
        result = {key: row[key] for key in keep if key in row}
        if row.get("operation") == "exact_positive_integer_expression":
            result.update(
                {
                    "opaque_definition": row.get("definition"),
                    "materialized": False,
                }
            )
        return result

    def _source_ref(self, ref: Any, row: Mapping[str, Any] | None = None) -> dict[str, Any]:
        node = self._ref_id(ref)
        if node is not None:
            item: dict[str, Any] = {
                "kind": "native_source_reference",
                "node": node,
                "materialized": False,
            }
            if node < self.prefix_upper:
                item["prefix_source_function"] = True
                self.opaque_inputs.add(str(node))
                self.counts["prefix_source_references"] += 1
            if node == self.N_node:
                item.update(
                    {
                        "name": "N",
                        "opaque_positive_integer": True,
                        "materialized": False,
                    }
                )
            if row is not None:
                operation = row.get("operation", row.get("op"))
                if operation is not None:
                    item["source_operation"] = operation
                provenance = self._provenance(row)
                if provenance:
                    item["provenance"] = provenance
            return item
        if _is_mapping(ref):
            result = {
                "kind": "opaque_source_expression_reference",
                "materialized": False,
                "provenance": self._provenance(ref),
            }
            self.opaque_inputs.add(self._ref_key(ref))
            return result
        return {
            "kind": "opaque_source_literal_reference",
            "literal": _scalar(ref),
            "materialized": False,
        }

    def _field_ref(self, raw: Any, *, literal_ok: bool = True) -> dict[str, Any]:
        if _is_mapping(raw):
            if any(key in raw for key in ("node", "node_id", "id", "source_node", "ref")):
                return self._compile_value(raw)
            if "operation" in raw or "op" in raw:
                return self._compile_inline_value(raw)
            return self._source_ref(raw)
        if raw is None:
            return self._source_ref(raw)
        if self._raw_row(raw) is not None or self._ref_key(raw) in self.record_keys:
            return self._compile_value(raw)
        if literal_ok:
            return self._source_ref(raw)
        raise KernelInputError("unresolved graph reference " + repr(raw))

    def _args(self, row: Mapping[str, Any]) -> list[Any]:
        args = row.get("arguments", row.get("children", row.get("terms")))
        if not isinstance(args, Sequence) or isinstance(args, (str, bytes)):
            raise KernelInputError("sum/product row lacks an arguments sequence")
        return list(args)

    def _binary_fields(self, row: Mapping[str, Any]) -> tuple[Any, Any]:
        if "numerator" in row and "denominator" in row:
            return row["numerator"], row["denominator"]
        if "left" in row and "right" in row:
            return row["left"], row["right"]
        raise KernelInputError("positive quotient row lacks numerator/denominator")

    def _new_value(
        self,
        operation: str,
        arguments: Sequence[Any] = (),
        *,
        source_node: Any = None,
        attributes: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        self._value_counter += 1
        ref = {"kind": "value_expression_ref", "id": "v" + str(self._value_counter)}
        row: dict[str, Any] = {
            "id": ref["id"],
            "kind": "value_expression",
            "operation": operation,
            "arguments": list(arguments),
            "numeric_evaluation": False,
        }
        if source_node is not None:
            row["source_node"] = self._ref_id(source_node)
        if attributes:
            row["attributes"] = dict(attributes)
        self.value_ir.append(row)
        self.counts["value_ir_rows"] += 1
        if source_node is not None:
            self.touched_nodes.add(self._ref_key(source_node))
        return ref

    def _new_tangent(
        self,
        operation: str,
        arguments: Sequence[Any] = (),
        *,
        source_node: Any = None,
        attributes: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        self._tangent_counter += 1
        ref = {"kind": "tangent_ref", "id": "dy" + str(self._tangent_counter)}
        row: dict[str, Any] = {
            "id": ref["id"],
            "kind": "tangent_expression",
            "operation": operation,
            "arguments": list(arguments),
            "operator": "formal_first_variation",
            "numeric_evaluation": False,
        }
        if source_node is not None:
            row["source_node"] = self._ref_id(source_node)
        if attributes:
            row["attributes"] = dict(attributes)
        self.tangent_ir.append(row)
        self.counts["tangent_ir_rows"] += 1
        return ref

    def _emit(self, event: str, row: Mapping[str, Any]) -> None:
        if event not in CALLBACK_EVENT_COUNTS:
            raise KernelInputError("callback event is outside the frozen audit contract: " + event)
        item = {"event": event, "row": dict(row)}
        self.events.append(item)
        self.counts["callback_events"] += 1
        if self.callback is not None:
            self.callback(event, item["row"])

    def _callback_event_counts(self) -> dict[str, int]:
        observed = {
            event: sum(1 for item in self.events if item.get("event") == event)
            for event in CALLBACK_EVENT_COUNTS
        }
        expected = dict(CALLBACK_EVENT_COUNTS)
        if observed != expected:
            raise KernelInputError(
                "callback audit event counts changed: expected "
                + repr(expected)
                + ", observed "
                + repr(observed)
            )
        return observed

    def _manifest_counts(self, event_counts: Mapping[str, int]) -> dict[str, int]:
        # Reconstruct the nonzero audit counters from the fixed callback rows,
        # then return the frozen literal contract.  This catches accidental
        # callback drift while keeping detailed compiler row counts in payload.
        derived = {name: 0 for name in SUCCESS_COUNTS}
        for event, targets in COUNT_EVENT_MAP.items():
            amount = int(event_counts.get(event, 0))
            for target in targets:
                if target not in derived:
                    raise KernelInputError("count event target is outside audit contract: " + str(target))
                derived[target] += amount
        if derived != dict(SUCCESS_COUNTS):
            raise KernelInputError(
                "callback audit counters changed: expected "
                + repr(SUCCESS_COUNTS)
                + ", derived "
                + repr(derived)
            )
        return dict(SUCCESS_COUNTS)

    def _equation(
        self,
        source_node: Any,
        operation: str,
        context: _Context,
        tangent: Mapping[str, Any],
        *,
        premise: str,
        extra: Mapping[str, Any] | None = None,
    ) -> None:
        row: dict[str, Any] = {
            "source_node": self._ref_id(source_node),
            "source_operation": operation,
            "mode": context.mode,
            "held_angle": context.held_angle,
            "dummy_variables": list(context.dummy_variables),
            "tangent_ref": tangent,
            "equation": premise,
            "evaluated": False,
        }
        if extra:
            row.update(extra)
        self.equations.append(row)
        self.counts["equation_receipts"] += 1

    def _compile_inline_value(self, row: Mapping[str, Any]) -> dict[str, Any]:
        operation = self._op(row, None)
        if operation == "exact_rational":
            return self._new_value(
                operation,
                (),
                attributes={
                    "numerator": _scalar(row.get("numerator")),
                    "denominator": _scalar(row.get("denominator")),
                },
            )
        raise KernelInputError("inline source operation is not a literal value: " + operation)

    def _compile_value(self, ref: Any) -> dict[str, Any]:
        key = self._ref_key(ref)
        cache_key = (key,)
        if cache_key in self.value_cache:
            return self.value_cache[cache_key]
        node = self._ref_id(ref)
        row = self._raw_row(ref)
        if node is None or row is None:
            result = self._source_ref(ref, row)
            self.value_cache[cache_key] = result
            return result
        self.touched_nodes.add(str(node))
        self.counts["source_graph_nodes_seen"] += 1

        # Every prefix node is opaque.  This preserves source-function status
        # and prevents a large upstream recipe from becoming local arithmetic.
        if node < self.prefix_upper:
            result = self._source_ref(node, row)
            self.value_cache[cache_key] = result
            return result

        if cache_key in self._value_active:
            raise KernelInputError("cycle in literal source graph at node " + str(node))
        self._value_active.add(cache_key)
        try:
            operation = self._op(row, node)
            if operation in {
                "exact_rational",
                "exact_positive_integer_expression",
                "current_original_source_parameter",
                "current_original_leading_function_recipe",
                "bound_variable",
            }:
                result = self._source_ref(node, row)
            elif operation in {"sum", "product"}:
                result = self._new_value(
                    operation,
                    [self._field_ref(item) for item in self._args(row)],
                    source_node=node,
                    attributes=self._provenance(row),
                )
            elif operation == "positive_quotient":
                numerator, denominator = self._binary_fields(row)
                result = self._new_value(
                    operation,
                    [self._field_ref(numerator), self._field_ref(denominator)],
                    source_node=node,
                    attributes=self._provenance(row),
                )
            elif operation == "negative":
                result = self._new_value(
                    operation,
                    [self._field_ref(row.get("argument"))],
                    source_node=node,
                    attributes=self._provenance(row),
                )
            elif operation == "analytic_unary":
                name = row.get("name")
                if name not in ANALYTIC_NAMES:
                    raise KernelInputError("unknown analytic unary " + repr(name))
                result = self._new_value(
                    operation,
                    [self._field_ref(row.get("argument"))],
                    source_node=node,
                    attributes={"name": name, **self._provenance(row)},
                )
            elif operation == "exact_original_source_flat_sigma":
                result = self._new_value(
                    operation,
                    [self._field_ref(row.get("argument"))],
                    source_node=node,
                    attributes=self._provenance(row),
                )
            elif operation == "original_lazy_flat_branch":
                result = self._new_value(
                    operation,
                    [
                        self._field_ref(row.get("active_body")),
                        self._field_ref(row.get("flat_value")),
                        self._field_ref(row.get("flat_predicate")),
                    ],
                    source_node=node,
                    attributes={
                        **self._provenance(row),
                        "branch_predicate_retained": True,
                        "active_body_not_evaluated_on_flat_branch": True,
                    },
                )
            elif operation == "definite_integral":
                result = self._new_value(
                    operation,
                    [
                        self._field_ref(row.get("integrand")),
                        self._field_ref(row.get("lower")),
                        self._field_ref(row.get("upper")),
                    ],
                    source_node=node,
                    attributes={
                        **self._provenance(row),
                        "variable": str(row.get("variable")),
                        "fixed_bound_leibniz_source": True,
                    },
                )
            elif operation == "substitute_original_inverse_angle":
                result = self._new_value(
                    operation,
                    [
                        self._field_ref(row.get("body")),
                        self._field_ref(row.get("inverse_angle")),
                    ],
                    source_node=node,
                    attributes={
                        **self._provenance(row),
                        "angle_bound_after_substitution": True,
                    },
                )
            elif operation == "original_monotone_phase_inverse":
                result = self._new_value(
                    operation,
                    [
                        self._field_ref(row.get("phase_function")),
                        self._field_ref(row.get("target_phase")),
                    ],
                    source_node=node,
                    attributes={
                        **self._provenance(row),
                        "angle_variable": str(row.get("angle_variable")),
                        "positive_inverse_derivative_premise": True,
                    },
                )
            elif operation == "function_substitution":
                result = self._new_value(
                    operation,
                    [
                        self._field_ref(row.get("expression")),
                        self._field_ref(row.get("variable")),
                        self._field_ref(row.get("value")),
                    ],
                    source_node=node,
                    attributes={
                        **self._provenance(row),
                        "Z_independent_substitution": bool(row.get("Z_independent_substitution")),
                    },
                )
            elif operation == "exact_real_comparison":
                result = self._new_value(
                    operation,
                    [self._field_ref(row.get("left")), self._field_ref(row.get("right"))],
                    source_node=node,
                    attributes={"operator": row.get("operator"), "predicate_not_differentiated": True},
                )
            else:
                raise KernelInputError("unsupported value operation " + operation)
            self.value_cache[cache_key] = result
            return result
        finally:
            self._value_active.discard(cache_key)

    def _is_angle_ref(self, ref: Any, name: str | None = None) -> bool:
        node = self._ref_id(ref)
        if node == self.angle_node:
            return True
        bound = self._bound_name(ref)
        return bound is not None and bound == (name or self.angle_name)

    def _is_exact_zero(self, ref: Any) -> bool:
        node = self._ref_id(ref)
        row = self._raw_row(ref)
        return bool(
            node is not None
            and row is not None
            and row.get("operation") == "exact_rational"
            and str(row.get("numerator")) == "0"
        )

    def _fixed_zero(self, reason: str, source_node: Any = None) -> dict[str, Any]:
        return self._new_tangent(
            "fixed_source_tangent",
            (),
            source_node=source_node,
            attributes={
                "equation": "D_y fixed source under declared local premise = 0",
                "reason": reason,
                "numeric_zero_not_materialized": True,
                "global_incoming_x_independence_proved": False,
            },
        )

    def _fixed_angle_zero(self, reason: str, source_node: Any = None) -> dict[str, Any]:
        return self._new_tangent(
            "typed_symbolic_zero",
            (),
            source_node=source_node,
            attributes={
                "equation": "angle partial of source terminal = 0",
                "reason": reason,
                "mode": "angle",
                "numeric_zero_not_materialized": True,
                "source_value_not_zeroed": True,
                "prefix_source_anglepartial_zero": True,
            },
        )

    def _one(self, source_node: Any, reason: str) -> dict[str, Any]:
        return self._new_tangent(
            "native_bound_variable_tangent",
            (),
            source_node=source_node,
            attributes={
                "equation": "D_y x=1" if reason == "native_x" else "D_angle angle=1",
                "reason": reason,
                "exact_symbolic_one": True,
            },
        )

    def _held_prefix_tangent(self, node: int, context: _Context) -> dict[str, Any]:
        if context.mode == "angle":
            return self._new_tangent(
                "typed_symbolic_zero",
                (),
                source_node=node,
                attributes={
                    "equation": "D_angle source_function_of_Z with local x held fixed = typed symbolic zero",
                    "source_domain_premise": True,
                    "value_retained_opaque": True,
                    "source_value_not_zeroed": True,
                    "local_x_held_fixed": True,
                    "global_incoming_x_independence_proved": False,
                    "prefix_source_anglepartial_zero": True,
                    "numeric_evaluation": False,
                },
            )
        return self._new_tangent(
            "typed_symbolic_zero",
            (),
            source_node=node,
            attributes={
                "equation": "D_y source_function_of_Z with local x held fixed = typed symbolic zero",
                "source_domain_premise": True,
                "value_retained_opaque": True,
                "source_value_not_zeroed": True,
                "local_x_held_fixed": True,
                "global_incoming_x_independence_proved": False,
                "numeric_evaluation": False,
            },
        )

    def _saved_history_tangent(self, node: int, context: _Context) -> dict[str, Any] | None:
        if context.mode == "angle":
            return self._fixed_angle_zero("source history has no physical angle variable", node)
        item = self.history.get(str(node), self.history.get(node))
        if not _is_mapping(item):
            return None
        expression = item.get("source_expression")
        if not _is_mapping(expression):
            raise KernelInputError("source history node " + str(node) + " lacks source_expression")
        if expression.get("operator") != "D_y":
            raise KernelInputError("source history node " + str(node) + " is not a saved D_y expression")
        if expression.get("evaluated") is True:
            raise KernelInputError("evaluated source history data is not admitted")
        if str(node) not in self.source_history_refs:
            tangent = self._new_tangent(
                "saved_source_differential_reference",
                (),
                source_node=node,
                attributes={
                    "field": item.get("field"),
                    "ordinary_Z_order": item.get("ordinary_Z_order"),
                    "proof_kind": item.get("proof_kind"),
                    "saved_result_pointer": item.get("saved_result_pointer"),
                    "source_expression": expression,
                    "provenance": "accepted saved source proof expression",
                    "numeric_evaluation": False,
                },
            )
            self.source_history_refs[str(node)] = tangent
            self.counts["source_history_joins"] += 1
        return self.source_history_refs[str(node)]

    def _t_sum(self, args: Sequence[Any], *, source_node: Any = None, equation: str = "sum rule") -> dict[str, Any]:
        return self._new_tangent("sum", list(args), source_node=source_node, attributes={"equation": equation})

    def _t_product(
        self,
        left_value: Any,
        right_value: Any,
        left_tangent: Any,
        right_tangent: Any,
        *,
        source_node: Any = None,
    ) -> dict[str, Any]:
        first = self._new_tangent(
            "product",
            [left_tangent, right_value],
            attributes={"factor": "left_tangent_times_right_value"},
        )
        second = self._new_tangent(
            "product",
            [left_value, right_tangent],
            attributes={"factor": "left_value_times_right_tangent"},
        )
        return self._t_sum(
            [first, second],
            source_node=source_node,
            equation="D_y(u*v)=D_y(u)*v+u*D_y(v)",
        )

    def _t_negative(self, tangent: Any, *, source_node: Any = None) -> dict[str, Any]:
        return self._new_tangent(
            "negative", [tangent], source_node=source_node, attributes={"equation": "D_y(-u)=-D_y(u)"}
        )

    def _t_quotient(
        self,
        numerator_value: Any,
        denominator_value: Any,
        numerator_tangent: Any,
        denominator_tangent: Any,
        *,
        source_node: Any = None,
        certificate: Any = None,
    ) -> dict[str, Any]:
        first = self._new_tangent(
            "product",
            [numerator_tangent, denominator_value],
            attributes={"equation": "D_y numerator_times_denominator"},
        )
        second = self._new_tangent(
            "product",
            [numerator_value, denominator_tangent],
            attributes={"equation": "numerator_times_D_y denominator"},
        )
        top = self._t_sum([first, self._t_negative(second)], equation="quotient numerator")
        bottom = self._new_value(
            "product",
            [denominator_value, denominator_value],
            attributes={"equation": "denominator_squared", "positive_denominator": certificate},
        )
        return self._new_tangent(
            "positive_quotient",
            [top, bottom],
            source_node=source_node,
            attributes={
                "equation": "D_y(u/v)=(D_y(u)*v-u*D_y(v))/v^2",
                "positive_denominator_certificate": certificate,
            },
        )

    def _compose_tangent_at_inverse(
        self,
        tangent: Any,
        inverse: Any,
        angle_name: str,
        *,
        partial: str,
    ) -> dict[str, Any]:
        inverse_value = self._compile_value(inverse)
        return self._new_tangent(
            "substitute_original_inverse_angle",
            [tangent, inverse_value],
            attributes={
                "angle_variable": angle_name,
                "partial": partial,
                "angle_bound_after_substitution": True,
                "bound_angle_not_free": True,
                "numeric_evaluation": False,
            },
        )

    def _compose_value_at_inverse(self, body: Any, inverse: Any, angle_name: str) -> dict[str, Any]:
        body_value = self._compile_value(body)
        inverse_value = self._compile_value(inverse)
        return self._new_value(
            "substitute_original_inverse_angle",
            [body_value, inverse_value],
            attributes={
                "angle_variable": angle_name,
                "angle_bound_after_substitution": True,
                "bound_angle_not_free": True,
            },
        )

    def _integrand_at_bound(
        self, integrand: Any, upper: Any, variable: str, *, lower: bool = False
    ) -> dict[str, Any]:
        integrand_value = self._compile_value(integrand)
        bound_value = self._compile_value(upper)
        return self._new_value(
            "bound_variable_substitution",
            [integrand_value, bound_value],
            attributes={
                "bound_variable": variable,
                "endpoint": "lower" if lower else "upper",
                "same_bound_dummy_not_differentiated": True,
                "numeric_evaluation": False,
            },
        )

    def _differentiate_integral(
        self, node: int, row: Mapping[str, Any], context: _Context
    ) -> dict[str, Any]:
        integrand = row.get("integrand")
        lower = row.get("lower")
        upper = row.get("upper")
        variable = str(row.get("variable"))
        at_upper = self._integrand_at_bound(integrand, upper, variable)

        own_upper_angle = (
            context.mode == "angle"
            and context.held_angle is None
            and variable not in context.dummy_variables
            and self._is_angle_ref(upper, variable)
        )
        if own_upper_angle:
            tangent = self._new_tangent(
                "definite_integral_upper_angle",
                [at_upper],
                source_node=node,
                attributes={
                    "equation": "D_angle integral_lower^angle f(s)ds=f(angle)",
                    "fixed_bound_leibniz": True,
                    "no_derivative_under_same_bound_dummy": True,
                    "bound_variable": variable,
                    "upper_endpoint_tangent": "exact_symbolic_one",
                    "lower_endpoint_tangent": "fixed_source_angle_partial",
                },
            )
            self.counts["integral_leibniz_rows"] += 1
            return tangent

        body_context = context.with_dummy(variable)
        body_tangent = self._differentiate(integrand, context=body_context)
        lower_value = self._compile_value(lower)
        upper_value = self._compile_value(upper)
        interior = self._new_tangent(
            "definite_integral",
            [body_tangent, lower_value, upper_value],
            source_node=node,
            attributes={
                "equation": "fixed-bound integral of held-angle integrand tangent",
                "integrand_tangent": body_tangent,
                "lower_bound_value": lower_value,
                "upper_bound_value": upper_value,
                "bound_variable": variable,
                "fixed_bound_leibniz": True,
                "no_derivative_under_same_bound_dummy": True,
            },
        )
        lower_tangent = self._differentiate(lower, context=context)
        upper_tangent = self._differentiate(upper, context=context)
        derivative_label = "D_" + context.mode
        upper_boundary = self._new_tangent(
            "product",
            [at_upper, upper_tangent],
            attributes={
                "equation": "integrand(upper)*" + derivative_label + " upper",
                "endpoint": "upper",
                "variation_mode": context.mode,
            },
        )
        at_lower = self._integrand_at_bound(integrand, lower, variable, lower=True)
        lower_boundary = self._new_tangent(
            "product",
            [at_lower, lower_tangent],
            attributes={
                "equation": "integrand(lower)*" + derivative_label + " lower",
                "endpoint": "lower",
                "variation_mode": context.mode,
            },
        )
        tangent = self._t_sum(
            [interior, upper_boundary, self._t_negative(lower_boundary)],
            source_node=node,
            equation="Leibniz fixed-bound rule with endpoint terms in " + context.mode + " mode",
        )
        self.counts["integral_leibniz_rows"] += 1
        return tangent

    def _differentiate_sigma(
        self, node: int, row: Mapping[str, Any], context: _Context
    ) -> dict[str, Any]:
        argument = row.get("argument")
        argument_value = self._compile_value(argument)
        argument_tangent = self._differentiate(argument, context=context)
        tangent = self._new_tangent(
            "exact_original_source_flat_sigma_variation",
            [argument_value, argument_tangent],
            source_node=node,
            attributes={
                "equation": "piecewise D_y sigma(s): 0 on flat extensions; sigma(s)(1-sigma(s))*(2/s^3+2/(1-s)^3)*D_y s on 0<s<1",
                "fresh_scalar_exact_definition": True,
                "active_open_interval": "0<s<1",
                "reciprocal_poles_excluded": True,
                "piecewise_exact_source_variation": True,
                "endpoint_branches": "s<=0 or s>=1",
                "endpoint_tangent": "typed_symbolic_zero",
                "piecewise_branches": [
                    {"condition": "s<=0", "tangent": "typed_symbolic_zero"},
                    {"condition": "0<s<1", "tangent": "sigma(s)*(1-sigma(s))*(2/s^3+2/(1-s)^3)*D_y s"},
                    {"condition": "s>=1", "tangent": "typed_symbolic_zero"},
                ],
                "dormant_reciprocal_poles": True,
                "no_zero_times_infinity": True,
                "reciprocal_expression_only_on_active_branch": True,
                "endpoint_formula_not_instantiated": True,
                "dummy_zero_not_substituted_into_reciprocals": True,
                "sigma_definition": row.get("definition"),
                "source_binding": row.get("source_binding"),
                "exact_function_not_its_interval_cap": row.get("exact_function_not_its_interval_cap"),
                "numeric_evaluation": False,
            },
        )
        self.counts["sigma_variation_rows"] += 1
        return tangent

    def _inverse_tangent(
        self, node: int, row: Mapping[str, Any], context: _Context
    ) -> dict[str, Any]:
        phase = row.get("phase_function")
        target = row.get("target_phase")
        angle_name = str(row.get("angle_variable"))
        target_tangent = self._differentiate(target, context=context)
        held_phase_tangent = self._differentiate(
            phase, context=context.with_held_angle(angle_name)
        )
        held_at_inverse = self._compose_tangent_at_inverse(
            held_phase_tangent,
            node,
            angle_name,
            partial="held_angle_D_y_Phi",
        )
        numerator = self._t_sum(
            [target_tangent, self._t_negative(held_at_inverse)],
            source_node=node,
            equation="D_"
            + context.mode
            + " target-[held-angle D_"
            + context.mode
            + " Phi](inverse)",
        )

        supplied_angle_derivative = None
        for key in ("phi_angle", "phase_angle_derivative", "angle_derivative"):
            if key in row:
                supplied_angle_derivative = row[key]
                break
        if supplied_angle_derivative is not None:
            denominator = self._compose_value_at_inverse(
                supplied_angle_derivative, node, angle_name
            )
            denominator_origin = "exact supplied Phi_angle source node"
        else:
            angle_tangent = self._differentiate(
                phase, context=_Context("angle", None, ())
            )
            denominator = self._compose_tangent_at_inverse(
                angle_tangent,
                node,
                angle_name,
                partial="Phi_angle",
            )
            denominator_origin = "fresh D_angle Phi source variation"

        certificate = row.get("derivative_positive_certificate")
        if not certificate:
            raise KernelInputError("monotone inverse lacks positive Phi_angle source premise")
        tangent = self._new_tangent(
            "positive_quotient",
            [numerator, denominator],
            source_node=node,
            attributes={
                "equation": "D_"
                + context.mode
                + " inverse=(D_"
                + context.mode
                + " target-[held-angle D_"
                + context.mode
                + " Phi](inverse))/Phi_angle(inverse)",
                "denominator": denominator,
                "Phi_angle_origin": denominator_origin,
                "positive_denominator_certificate": certificate,
                "positive_inverse_derivative_conditional": True,
                "N_materialized": False,
                "numeric_evaluation": False,
            },
        )
        self.counts["inverse_ift_rows"] += 1
        self.flags["conditional_positive_inverse_derivative"] = True
        return tangent

    def _substitute_tangent(
        self, node: int, row: Mapping[str, Any], context: _Context
    ) -> dict[str, Any]:
        body = row.get("body")
        inverse = row.get("inverse_angle")
        angle_name = str(row.get("angle_variable"))
        direct = self._differentiate(
            body,
            context=context.with_held_angle(angle_name),
        )
        direct_at_inverse = self._compose_tangent_at_inverse(
            direct,
            inverse,
            angle_name,
            partial="direct_partial_at_inverse",
        )
        body_angle = self._differentiate(
            body,
            context=_Context("angle", None, context.dummy_variables),
        )
        body_angle_at_inverse = self._compose_tangent_at_inverse(
            body_angle,
            inverse,
            angle_name,
            partial="body_angle",
        )
        inverse_tangent = self._differentiate(inverse, context=context)
        chain = self._new_tangent(
            "product",
            [body_angle_at_inverse, inverse_tangent],
            attributes={
                "equation": "body_angle(inverse)*D_" + context.mode + " inverse",
                "bound_angle_not_free": True,
                "variation_mode": context.mode,
            },
        )
        tangent = self._t_sum(
            [direct_at_inverse, chain],
            source_node=node,
            equation="substitute inverse chain rule: direct partial+body_angle*inverse tangent",
        )
        return tangent

    def _differentiate(self, ref: Any, *, context: _Context) -> dict[str, Any]:
        key = (self._ref_key(ref),) + context.key()
        if key in self.tangent_cache:
            return self.tangent_cache[key]
        node = self._ref_id(ref)
        row = self._raw_row(ref)
        if node is None or row is None:
            if node is not None and node < self.prefix_upper:
                result = self._held_prefix_tangent(node, context)
            else:
                raise KernelInputError("cannot differentiate opaque non-prefix reference " + repr(ref))
            self.tangent_cache[key] = result
            return result

        # Accepted m/e/p ordinary source rows carry prior source-proof
        # differentials.  They are DATA references, not a recurrence to replay.
        if str(node) in {str(item) for item in self.history}:
            result = self._saved_history_tangent(node, context)
            if result is None:
                raise KernelInputError("source history tangent missing for node " + str(node))
            self.tangent_cache[key] = result
            self._equation(
                node,
                str(row.get("operation", "source_history")),
                context,
                result,
                premise="saved accepted source proof expression",
            )
            return result

        # A graph row below the local chart boundary is a source function,
        # never an arithmetic recipe for this compiler.
        if node < self.prefix_upper:
            saved = self._saved_history_tangent(node, context)
            if saved is not None:
                result = saved
            elif node == self.N_node:
                result = (
                    self._fixed_angle_zero("N has no angle dependence", node)
                    if context.mode == "angle"
                    else self._fixed_zero("N is a fixed opaque source parameter", node)
                )
            elif row.get("operation") in {
                "exact_rational",
                "exact_positive_integer_expression",
                "current_original_source_parameter",
            }:
                result = (
                    self._fixed_angle_zero("fixed source parameter", node)
                    if context.mode == "angle"
                    else self._fixed_zero("fixed source parameter", node)
                )
            else:
                result = self._held_prefix_tangent(node, context)
            self.tangent_cache[key] = result
            self._equation(
                node,
                str(row.get("operation", "opaque_source")),
                context,
                result,
                premise="prefix source retained under local source-domain premise",
            )
            return result

        if key in self._tangent_active:
            raise KernelInputError("cycle in tangent graph at node " + str(node))
        self._tangent_active.add(key)
        try:
            operation = self._op(row, node)
            if operation in {
                "exact_rational",
                "exact_positive_integer_expression",
                "current_original_source_parameter",
            }:
                result = (
                    self._fixed_angle_zero("fixed source parameter", node)
                    if context.mode == "angle"
                    else self._fixed_zero("fixed source parameter", node)
                )
                premise = "fixed source parameter tangent"
            elif operation == "bound_variable":
                name = str(row.get("name"))
                if name in context.dummy_variables:
                    result = self._fixed_zero("bound integral dummy held fixed", node)
                    premise = "bound dummy variable held fixed under Leibniz rule"
                elif node == self.x_node:
                    result = self._one(node, "native_x") if context.mode == "y" else self._fixed_angle_zero("x has no angle dependence", node)
                    premise = "native y=D_x"
                elif node == self.angle_node or name == self.angle_name:
                    if context.mode == "angle" and context.held_angle != name:
                        result = self._one(node, "angle")
                    elif context.mode == "angle":
                        result = self._fixed_angle_zero("angle held in IFT partial", node)
                    else:
                        result = self._fixed_zero("angle held in total-y partial", node)
                    premise = "free angle variable / held-angle distinction"
                elif name == "Z":
                    result = self._held_prefix_tangent(node, context)
                    premise = "ordinary Z source variable held in local radial tangent"
                else:
                    raise KernelInputError("unclassified bound variable " + name)
            elif operation == "sum":
                args = self._args(row)
                result = self._t_sum(
                    [self._differentiate(arg, context=context) for arg in args],
                    source_node=node,
                    equation="sum rule",
                )
                premise = "sum rule"
            elif operation == "product":
                args = self._args(row)
                if len(args) < 2:
                    raise KernelInputError("product row has fewer than two factors at node " + str(node))
                result = self._differentiate(args[0], context=context)
                if len(args) == 2:
                    result = self._t_product(
                        self._compile_value(args[0]),
                        self._compile_value(args[1]),
                        result,
                        self._differentiate(args[1], context=context),
                        source_node=node,
                    )
                else:
                    current_value = self._compile_value(args[0])
                    current_tangent = result
                    for arg in args[1:]:
                        current_tangent = self._t_product(
                            current_value,
                            self._compile_value(arg),
                            current_tangent,
                            self._differentiate(arg, context=context),
                            source_node=node,
                        )
                        current_value = self._new_value(
                            "product",
                            [current_value, self._compile_value(arg)],
                            attributes={"folded_source_product": True},
                        )
                    result = current_tangent
                premise = "product rule"
            elif operation == "positive_quotient":
                numerator, denominator = self._binary_fields(row)
                result = self._t_quotient(
                    self._compile_value(numerator),
                    self._compile_value(denominator),
                    self._differentiate(numerator, context=context),
                    self._differentiate(denominator, context=context),
                    source_node=node,
                    certificate=row.get("source_positive_certificate"),
                )
                premise = "positive quotient rule"
            elif operation == "negative":
                result = self._t_negative(
                    self._differentiate(row.get("argument"), context=context),
                    source_node=node,
                )
                premise = "negative rule"
            elif operation == "analytic_unary":
                name = row.get("name")
                if name not in ANALYTIC_NAMES:
                    raise KernelInputError("unknown analytic unary " + repr(name))
                argument = row.get("argument")
                argument_value = self._compile_value(argument)
                unary_value = self._compile_value(node)
                argument_tangent = self._differentiate(argument, context=context)
                if name == "exp":
                    result = self._new_tangent(
                        "product",
                        [unary_value, argument_tangent],
                        source_node=node,
                        attributes={"equation": "D_y exp(u)=exp(u)*D_y u"},
                    )
                elif name == "log":
                    result = self._new_tangent(
                        "positive_quotient",
                        [argument_tangent, argument_value],
                        source_node=node,
                        attributes={
                            "equation": "D_y log(u)=D_y u/u",
                            "positive_denominator_certificate": "positive log argument source premise",
                        },
                    )
                elif name == "positive_sqrt":
                    two = self._compile_value(48)
                    denominator = self._new_value(
                        "product",
                        [two, unary_value],
                        attributes={"equation": "2*positive_sqrt_argument"},
                    )
                    result = self._new_tangent(
                        "positive_quotient",
                        [argument_tangent, denominator],
                        source_node=node,
                        attributes={
                            "equation": "D_y sqrt(u)=D_y u/(2*sqrt(u))",
                            "positive_denominator_certificate": "positive sqrt source premise",
                        },
                    )
                elif name in {"cos", "sin"}:
                    companion_value = self._new_value(
                        "analytic_unary",
                        [argument_value],
                        attributes={"name": "sin" if name == "cos" else "cos", "derived_companion": True},
                    )
                    if name == "cos":
                        negative_companion = self._new_value(
                            "negative",
                            [companion_value],
                            attributes={"equation": "-sin(argument)"},
                        )
                        result = self._new_tangent(
                            "product",
                            [negative_companion, argument_tangent],
                            source_node=node,
                            attributes={"equation": "D_y cos(u)=-sin(u)*D_y u"},
                        )
                    else:
                        result = self._new_tangent(
                            "product",
                            [companion_value, argument_tangent],
                            source_node=node,
                            attributes={"equation": "D_y sin(u)=cos(u)*D_y u"},
                        )
                else:
                    result = self._new_tangent(
                        "fractional_part_variation",
                        [argument_tangent],
                        source_node=node,
                        attributes={
                            "equation": "D_y frac(u)=D_y u on an open off-seam phase cell",
                            "off_phase_seam": True,
                            "phase_node": node == self.phase_node,
                            "N_materialized": False,
                        },
                    )
                    self.counts["phase_fractional_rows"] += 1
                    self.flags["actual_phase_fractional_part_tangent"] |= node == self.phase_node
                premise = "analytic unary rule " + str(name)
            elif operation == "exact_original_source_flat_sigma":
                result = self._differentiate_sigma(node, row, context)
                premise = "fresh exact sigma-definition variation on active open interval"
            elif operation == "original_lazy_flat_branch":
                active = row.get("active_body")
                flat = row.get("flat_value")
                predicate = self._compile_value(row.get("flat_predicate"))
                active_tangent = self._differentiate(active, context=context)
                if self._is_exact_zero(flat):
                    flat_tangent = self._fixed_zero("native flat_value is exact source zero", node)
                    flat_kind = "flat_value_native_zero"
                else:
                    flat_tangent = self._differentiate(flat, context=context)
                    flat_kind = "flat_value_variation"
                result = self._new_tangent(
                    "original_lazy_flat_branch",
                    [active_tangent, flat_tangent, predicate],
                    source_node=node,
                    attributes={
                        "equation": "branchwise tangent; predicate retained",
                        "branch_predicate": predicate,
                        "flat_tangent_kind": flat_kind,
                        "active_branch_variation": active_tangent,
                        "flat_branch_variation": flat_tangent,
                        "active_body_not_evaluated_on_flat_branch": True,
                        "no_global_seam_claim": True,
                    },
                )
                self.counts["lazy_branch_rows"] += 1
                premise = "open active/flat branch with predicate preserved"
            elif operation == "definite_integral":
                result = self._differentiate_integral(node, row, context)
                premise = "fixed-bound Leibniz rule"
            elif operation == "substitute_original_inverse_angle":
                result = self._substitute_tangent(node, row, context)
                premise = "inverse-angle chain rule"
            elif operation == "original_monotone_phase_inverse":
                result = self._inverse_tangent(node, row, context)
                premise = "implicit-function tangent with positive Phi_angle"
            elif operation == "function_substitution":
                raise KernelInputError(
                    "local function_substitution tangent is outside the bounded source contract"
                )
            elif operation == "exact_real_comparison":
                raise KernelInputError("predicate node " + str(node) + " was requested for differentiation")
            elif operation == "current_original_leading_function_recipe":
                result = self._held_prefix_tangent(node, context)
                premise = "opaque source recipe tangent"
            else:
                raise KernelInputError("unhandled tangent operation " + operation)

            self.tangent_cache[key] = result
            self._equation(node, operation, context, result, premise=premise)
            return result
        finally:
            self._tangent_active.discard(key)

    def _root(self, name: str) -> Any:
        if name in self.roots:
            return self.roots[name]
        raise KernelInputError("required root " + name + " is missing")

    def _binding(self, section: str, name: str) -> Any:
        bindings = self.capsule.get("source_bindings", {})
        if _is_mapping(bindings) and _is_mapping(bindings.get(section)):
            section_map = bindings[section]
            if name in section_map:
                return section_map[name]
        raise KernelInputError("source binding " + section + "." + name + " is missing")

    def _target_row(self, name: str, ref: Any, context: _Context) -> dict[str, Any]:
        value = self._compile_value(ref)
        tangent = self._differentiate(ref, context=context)
        raw = self._raw_row(ref)
        row = {
            "name": name,
            "source_node": self._ref_id(ref),
            "value_ref": value,
            "tangent_ref": tangent,
            "operator": "D_y",
            "ordinary_Z_order": raw.get("ordinary_Z_order") if raw else None,
            "evaluated": False,
            "numeric_evaluation": False,
        }
        return row

    def _source_history_rows(self) -> list[dict[str, Any]]:
        result = []
        keys = sorted(
            self.history,
            key=lambda item: int(item) if str(item).lstrip("-").isdigit() else str(item),
        )
        for raw_node in keys:
            node = _node_number(raw_node)
            if node is None:
                raise KernelInputError("source history key is not a native node: " + str(raw_node))
            tangent = self._saved_history_tangent(node, _Context("y"))
            if tangent is None:
                raise KernelInputError("source history tangent missing for node " + str(node))
            item = self.history[raw_node]
            result.append(
                {
                    "node": node,
                    "field": item.get("field"),
                    "ordinary_Z_order": item.get("ordinary_Z_order"),
                    "saved_result_pointer": item.get("saved_result_pointer"),
                    "tangent_ref": tangent,
                    "source_expression": item.get("source_expression"),
                    "angle_partial": "fixed_source_angle_partial",
                    "evaluated": False,
                }
            )
        return result

    def run(self) -> dict[str, Any]:
        source_history_rows = self._source_history_rows()
        y_context = _Context("y")

        target_refs: list[tuple[str, Any]] = [
            ("E_background", self._root("E_background")),
            ("V_background", self._root("V_background")),
            ("J_sigma", 241546),
            ("actual_phase", self._root("actual_phase")),
            ("Phi", self._root("Phi")),
            ("inverse_angle", self._binding("loop", "inverse")),
            ("phase_lambda", self._root("phase_lambda")),
            ("phase_lambda_at_inverse", self._root("phase_lambda_at_inverse")),
            ("A", self._root("A")),
            ("B_over_Pstar", self._root("B_over_Pstar")),
            ("E_N", self._root("E_N")),
            ("V_N", self._root("V_N")),
        ]
        targets = [self._target_row(name, ref, y_context) for name, ref in target_refs]

        self._emit(
            "source_history_join_summary",
            {
                "count": len(source_history_rows),
                "fields": sorted({str(item.get("field")) for item in source_history_rows}),
                "ordinary_Z_orders": sorted(
                    {
                        item.get("ordinary_Z_order")
                        for item in source_history_rows
                        if item.get("ordinary_Z_order") is not None
                    }
                ),
                "evaluated": False,
            },
        )
        phase_target = next(item for item in targets if item["name"] == "actual_phase")
        self._emit(
            "phase_variation_summary",
            {
                "source_node": phase_target["source_node"],
                "tangent_ref": phase_target["tangent_ref"],
                "fractional_part_rows": self.counts["phase_fractional_rows"],
                "off_phase_seam": True,
                "actual_affine_graph_expanded": True,
                "N_materialized": False,
                "evaluated": False,
            },
        )
        self._emit(
            "branch_integral_sigma_summary",
            {
                "integral_leibniz_rows": self.counts["integral_leibniz_rows"],
                "sigma_variation_rows": self.counts["sigma_variation_rows"],
                "lazy_branch_rows": self.counts["lazy_branch_rows"],
                "dummy_capture": True,
                "predicate_preserved": True,
                "evaluated": False,
            },
        )
        self._emit(
            "inverse_ift_summary",
            {
                "inverse_ift_rows": self.counts["inverse_ift_rows"],
                "fresh_Phi_angle_variation": True,
                "saved_phase_lambda_used_as_positive_provenance_only": True,
                "old_lambda_equality_claimed": False,
                "evaluated": False,
            },
        )

        # Retain the finite-N definitions as equations over the actual source
        # roots.  The tangent rows above come from recursively expanding those
        # roots; no stored Ey/Vy alias is consumed.
        N_value = self._source_ref(self.N_node, self._raw_row(self.N_node))
        E_value = self._compile_value(self._root("E_background"))
        V_value = self._compile_value(self._root("V_background"))
        A_value = self._compile_value(self._root("A"))
        B_value = self._compile_value(self._root("B_over_Pstar"))
        finite_definition = {
            "E_N": {
                "equation": "E_N=E_background*exp(A/N)",
                "actual_root": self._ref_id(self._root("E_N")),
                "E": E_value,
                "A": A_value,
                "N": N_value,
                "N_materialized": False,
            },
            "V_N": {
                "equation": "V_N=V_background+B_over_Pstar/N",
                "actual_root": self._ref_id(self._root("V_N")),
                "V": V_value,
                "B_over_Pstar": B_value,
                "N": N_value,
                "physical_B_not_substituted": True,
                "N_materialized": False,
            },
        }
        self._emit("finite_N_definition", finite_definition)

        self.flags["actual_E_background_exp_J_sigma_tangent"] = any(
            item["name"] in {"E_background", "J_sigma"} for item in targets
        )
        self.flags["actual_A_tangent"] = True
        self.flags["actual_B_over_Pstar_tangent"] = True
        self.flags["actual_E_N_tangent"] = True
        self.flags["actual_V_N_tangent"] = True
        self.flags["conditional_C2_regular_open_cell"] = True
        self.flags["conditional_native_O3_first_variation_compiled"] = True
        self._emit(
            "family_completed",
            {
                "family": FAMILY,
                "source_only": True,
                "numeric_evaluation": False,
                "N_materialized": False,
                "conditional_scope": "C2 regular open O3 phase cells and open lazy branches",
                "counts": dict(self.counts),
                "flags": dict(self.flags),
            },
        )
        self.counts["callback_events"] = len(self.events)
        callback_event_counts = self._callback_event_counts()
        actual_counts = self._manifest_counts(callback_event_counts)
        diagnostic_flags = dict(self.flags)
        return {
            "schema": SCHEMA,
            "component": COMPONENT,
            "family": FAMILY,
            "attempt_status": "passed",
            "passed": True,
            "metadata_only": False,
            "conditional": True,
            "conditional_on_declared_source_atoms": True,
            "source_only": True,
            "full_goal_complete": False,
            "numeric_evaluation": False,
            "runtime_model_metadata": None,
            "actual_counts": actual_counts,
            "new_flags": dict(SUCCESS_FLAGS),
            "proof_count": 1,
            "callback_event_count": sum(callback_event_counts.values()),
            "callback_event_counts": callback_event_counts,
            "count_event_map": {
                name: list(targets) for name, targets in COUNT_EVENT_MAP.items()
            },
            "flags": dict(SUCCESS_FLAGS),
            "failure": None,
            "domain_premises": {
                "regularity": "conditional C2 regular open O3 cells only",
                "x": "0<x<1",
                "phase": "actual fractional_part phase off seams",
                "y": "native y=log R and D_y=D_x",
                "inverse": "Phi_angle>0 from actual source recipe; no global conditioning claim",
                "sigma": "active 0<s<1 branch; reciprocal poles excluded",
                "lazy": "branchwise predicate-preserving variation",
                "incoming": "prefix nodes are opaque source functions of Z under local x-held-fixed premise",
                "incoming_global_independence": False,
                "theta_phi_distinct": True,
            },
            "payload": {
                "targets": targets,
                "finite_N_definition": finite_definition,
                "source_history_tangents": source_history_rows,
                "saved_flux_result_contract": {
                    "provided_as_mapping": bool(self.saved),
                    "schema": self.saved.get("schema") if _is_mapping(self.saved) else None,
                    "source_only": self.saved.get("source_only") if _is_mapping(self.saved) else None,
                    "numeric_evaluation": self.saved.get("numeric_evaluation") if _is_mapping(self.saved) else None,
                    "accepted_source_expression_data_referenced": True,
                    "old_runtime_not_replayed": True,
                },
                "source_premises": list(self.premises),
                "diagnostic_flags": diagnostic_flags,
                "internal_counts": dict(self.counts),
                "success_counts_contract": dict(SUCCESS_COUNTS),
                "success_flags_contract": dict(SUCCESS_FLAGS),
                "callback_event_counts_contract": dict(callback_event_counts),
                "count_event_map_contract": {
                    name: list(targets) for name, targets in COUNT_EVENT_MAP.items()
                },
                "value_ir": list(self.value_ir),
                "tangent_ir": list(self.tangent_ir),
                "equation_receipts": list(self.equations),
                "opaque_source_inputs": sorted(
                    self.opaque_inputs,
                    key=lambda item: int(item) if item.lstrip("-").isdigit() else item,
                ),
                "source_graph_nodes_touched": sorted(
                    self.touched_nodes,
                    key=lambda item: int(item) if item.lstrip("-").isdigit() else item,
                ),
                "saved_phase_lambda_provenance": {
                    "source_node": self._root("phase_lambda"),
                    "at_inverse_source_node": self._root("phase_lambda_at_inverse"),
                    "used_as_positive_source_premise_only": True,
                    "fresh_Phi_angle_variation_used_for_IFT": True,
                    "old_lambda_equality_claimed": False,
                },
                "events": list(self.events),
            },
        }


def _rejected_report(
    capsule: Any,
    error: Exception,
    callback: Callable[[str, Mapping[str, Any]], Any] | None,
) -> dict[str, Any]:
    failure = {
        "kind": "kernel_input_rejection",
        "error_type": type(error).__name__,
        "error": str(error),
        "callback_prefix": [],
        "stopped": True,
    }
    row = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "family": FAMILY,
        "attempt_status": "rejected",
        "passed": False,
        "metadata_only": False,
        "conditional": False,
        "source_only": True,
        "full_goal_complete": False,
        "numeric_evaluation": False,
        "runtime_model_metadata": None,
        "proof_count": 0,
        "failure": failure,
        "flags": {name: False for name in SUCCESS_FLAGS},
        "new_flags": {name: False for name in SUCCESS_FLAGS},
        "counts": {
            "unknown_operation_refusals": 1 if "unknown source operation" in str(error) else 0
        },
        "payload": {
            "events": [],
            "source_capsule_schema": capsule.get("schema") if _is_mapping(capsule) else None,
            "failure": failure,
        },
    }
    return row


def derive(
    capsule: Mapping[str, Any],
    saved_flux_result: Mapping[str, Any] | None,
    callback: Callable[[str, Mapping[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Compile first total-y variations from the supplied literal source graph.

    The function is source-only: all values and tangents are references or IR
    equations, ``N`` remains opaque, and no source runtime or numerical oracle
    is consulted.  Unknown operations and missing non-prefix source inputs are
    rejected in the returned report.
    """

    try:
        return _Compiler(capsule, saved_flux_result, callback).run()
    except KernelInputError as error:
        return _rejected_report(capsule, error, callback)


__all__ = [
    "FAMILY",
    "SCHEMA",
    "COMPONENT",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "CALLBACK_EVENT_COUNTS",
    "COUNT_EVENT_MAP",
    "KernelInputError",
    "derive",
]
