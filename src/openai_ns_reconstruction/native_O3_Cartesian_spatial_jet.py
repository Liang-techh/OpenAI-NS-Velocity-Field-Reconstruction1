"""Fresh source-only Cartesian spatial jet compiler for the O3 branch.

This module consumes literal saved-result data and declarative source-partial
recipes.  It builds a new scalar IR and runs a reverse adjoint pass over only
the saved, uncomposed Cartesian roots.  It never imports a source module,
re-enters an old differentiator, materializes ``N``, or evaluates a source
value.

The root runner owns freezing, entry, invocation, and saving.  Workers only
provide this unexecuted compiler and static evidence about its source.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "native_O3_Cartesian_spatial_jet_kernel_result_v1"
COMPONENT = "native_O3_Cartesian_spatial_jet"


# Literal contracts used by the root manifest extractor.  Per-branch counts
# are diagnostic payload data; this bounded contract counts method families.
SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "api_refusal": 0,
    "Cartesian_composition": 1,
    "curl": 1,
    "divergence": 1,
    "failed_family_calls": 0,
    "frame_derivation": 1,
    "geometry_dual_calls": 0,
    "new_algebra_family_calls": 1,
    "new_differential_family_calls": 1,
    "new_mathematical_proofs": 1,
    "numeric_evaluation_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "physical_gradient": 1,
    "saved_geometry_differentials_consumed": 1,
    "provider_calls": 0,
    "reverse_adjoint": 1,
    "source_family_calls": 1,
    "source_join": 1,
    "source_value_reducer_calls": 0,
    "temporal_recursion_calls": 0,
    "volume_integrand": 1,
}


SUCCESS_FLAGS = {
    "Cartesian_spatial_gradients_defined": True,
    "curl_defined": True,
    "source_volume_integrands_defined": True,
    "physical_volume_speed_integrand_defined": True,
    "physical_volume_vorticity_integrand_defined": True,
    "conditional_cartesian_spatial_gradients_defined": True,
    "conditional_curl_defined": True,
    "conditional_divergence_expression_defined": True,
    "conditional_physical_volume_integrands_defined": True,
    "conditional_native_x_volume_integrands_defined": True,
    "conditional_source_reverse_adjoint_defined": True,
    "conditional_saved_geometry_differentials_consumed": True,
    "conditional_saved_positive_inverse_composition_defined": True,
    "conditional_native_Rm_Ra_phase_Z_invariance_proved": True,
    "conditional_fixed_x_Z_equals_fixed_R_Z": True,
    "actual_source_oracles_available": False,
    "source_oracles": False,
    "actual_source_Z_derivative_admission": False,
    "source_frame_admitted_unconditionally": False,
    "actual_incoming_x_independence_proved": False,
    "global_incoming_x_independence_proved": False,
    "global_phase_seam_regularity_proved": False,
    "global_axis_regularity_proved": False,
    "axis_extension_proved": False,
    "seam_regularity_proved": False,
    "global_independence": False,
    "unconditional_Z4_source_admission": False,
    "Z4_admission": False,
    "Rm_Z_certificate_installed": False,
    "unconditional_physical_identity_proved": False,
    "unconditional_physical": False,
    "NS_identity_proved": False,
    "NS_proved": False,
    "temporal_recursion": False,
    "temporal_recursion_proved": False,
    "cap_deletion": False,
    "caps_proved": False,
    "numeric_evaluation": False,
    "N_materialized": False,
    "source_value_identity": False,
    "full_goal_complete": False,
}


# Exactly nine events are emitted on success.  Composition is recorded in the
# physical_gradients row, keeping the callback surface fixed and auditable.
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "source_frame_zero_differential": 1,
    "source_join": 1,
    "reverse_adjoint": 1,
    "geometry_data": 1,
    "physical_gradients": 1,
    "curl_divergence": 1,
    "volume_integrands": 1,
    "family_completed": 1,
}


COUNT_EVENT_MAP = {
    "family_entered": [
        "new_differential_family_calls",
        "new_algebra_family_calls",
        "source_family_calls",
    ],
    "source_frame_zero_differential": ["frame_derivation"],
    "source_join": ["source_join"],
    "reverse_adjoint": ["reverse_adjoint"],
    "geometry_data": ["saved_geometry_differentials_consumed"],
    "physical_gradients": ["physical_gradient", "Cartesian_composition"],
    "curl_divergence": ["curl", "divergence"],
    "volume_integrands": ["volume_integrand"],
    "family_completed": ["new_mathematical_proofs"],
}


SUPPORTED_REACHABLE_OPERATIONS = frozenset(
    {
        "rational",
        "native_function_reference",
        "free_variable",
        "bound_variable",
        "add",
        "multiply",
        "negative",
        "power",
        "quotient",
        "log",
        "exp",
        "cos",
        "sin",
    }
)
ANALYTIC_OPERATIONS = frozenset({"log", "exp", "cos", "sin"})
AXES = ("R", "Z", "theta")
VELOCITY_FIELDS = ("u", "v", "w")


# Legacy table mode may use only these exact fixed-port defaults.  Recipe mode
# below disables these defaults: every native source node then needs an
# explicit declarative recipe.  Rm is deliberately absent in both modes.
FIXED_PORT_PARTIALS = {
    "241787": {"R": 1, "Z": 0, "theta": 0},
    "20": {"R": 0, "Z": 1, "theta": 0},
    "33572": {"R": 0, "Z": 0, "theta": 0},
    "25220": {"R": 0, "Z": 0, "theta": 0},
}
FRAME_NAME = "fixed_profile_R_under_admitted_frame_contract"


class KernelInputError(ValueError):
    """Raised when a literal input violates this narrow source-only contract."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _saved_id(value: Any) -> str | None:
    if isinstance(value, Mapping):
        for key in ("saved_expression_ref", "expression_ref", "saved_pointer", "ref", "id"):
            if key in value:
                found = _saved_id(value[key])
                if found is not None:
                    return found
        return None
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("physical_expr:"):
            return text
        if text.lstrip("-").isdigit():
            return "physical_expr:" + str(int(text))
    if isinstance(value, int) and not isinstance(value, bool):
        return "physical_expr:" + str(value)
    return None


def _native_id(value: Any) -> str:
    if isinstance(value, bool):
        raise KernelInputError("boolean is not a native source id")
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str) and value.strip().lstrip("-").isdigit():
        return str(int(value.strip()))
    raise KernelInputError("native_node must be an integer-like id")


def _json_scalar(value: Any) -> Any:
    return value if value is None or isinstance(value, (str, int, float, bool)) else str(value)


def _exact_fraction(value: Any) -> Fraction | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str) and value.strip().lstrip("-").isdigit():
        return Fraction(int(value.strip()), 1)
    if _is_mapping(value) and value.get("operation") in {"rational", "exact_rational"}:
        numerator, denominator = value.get("numerator"), value.get("denominator", 1)
        if isinstance(numerator, int) and isinstance(denominator, int) and denominator:
            return Fraction(numerator, denominator)
    return None


def _compact(value: Any) -> str:
    return str(value).replace("_", "").replace(" ", "").lower()


class _IRBuilder:
    """Translate saved value rows into fresh scalar IR rows."""

    def __init__(self, rows: Sequence[Mapping[str, Any]]) -> None:
        self.saved_rows = list(rows)
        self.saved_by_id: dict[str, Mapping[str, Any]] = {}
        for row in self.saved_rows:
            item = _as_mapping(row, "saved expression row")
            row_id = item.get("id")
            if not isinstance(row_id, str) or not row_id.startswith("physical_expr:"):
                raise KernelInputError("saved expression row id must be physical_expr:<n>")
            if row_id in self.saved_by_id:
                raise KernelInputError("duplicate saved expression id " + row_id)
            self.saved_by_id[row_id] = item
        self.saved_order = {row_id: i for i, row_id in enumerate(self.saved_by_id)}
        self.rows: list[dict[str, Any]] = []
        self.cache: dict[str, str] = {}
        self.primitive_cache: dict[str, str] = {}
        self.coordinate_cache: dict[str, str] = {}

    def _args(self, row: Mapping[str, Any]) -> list[str]:
        raw = row.get("arguments")
        if raw is None and "argument" in row:
            raw = [row["argument"]]
        if raw is None:
            return []
        if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray)):
            raise KernelInputError("saved arguments must be a sequence")
        result: list[str] = []
        for arg in raw:
            arg_id = _saved_id(arg)
            if arg_id is None:
                raise KernelInputError("saved argument is not a physical expression reference")
            result.append(arg_id)
        return result

    def row(self, ref: Any) -> Mapping[str, Any]:
        row_id = _saved_id(ref)
        if row_id is None or row_id not in self.saved_by_id:
            raise KernelInputError("unknown saved expression reference " + repr(ref))
        return self.saved_by_id[row_id]

    def row_id(self, ref: Any) -> str:
        row_id = _saved_id(ref)
        if row_id is None or row_id not in self.saved_by_id:
            raise KernelInputError("unknown saved expression reference " + repr(ref))
        return row_id

    def new(self, operation: str, arguments: Sequence[Any] = (), attrs: Mapping[str, Any] | None = None) -> str:
        row_id = "spatial_expr:" + str(len(self.rows))
        row: dict[str, Any] = {
            "id": row_id,
            "operation": operation,
            "arguments": list(arguments),
            "source_values_unevaluated": True,
            "numeric_evaluation": False,
        }
        if attrs:
            row.update(dict(attrs))
        self.rows.append(row)
        return row_id

    def rational(self, numerator: int, denominator: int = 1, *, premise: str | None = None) -> str:
        if not isinstance(numerator, int) or not isinstance(denominator, int) or denominator == 0:
            raise KernelInputError("malformed fresh rational")
        attrs: dict[str, Any] = {"numerator": numerator, "denominator": denominator}
        if premise is not None:
            attrs["premise"] = premise
        return self.new("rational", (), attrs)

    def _primitive(self, row_id: str, row: Mapping[str, Any]) -> str:
        if row_id in self.primitive_cache:
            return self.primitive_cache[row_id]
        native = _native_id(row.get("native_node"))
        path = row.get("path")
        text = _compact(path)
        for forbidden in ("p1", "m15", "p2mix"):
            if forbidden in text:
                raise KernelInputError("forbidden source lineage token " + forbidden)
        native_value: Any = int(native) if native.lstrip("-").isdigit() else native
        result = self.new(
            "primitive_source_port",
            (),
            {
                "opaque_id": "native:" + native + ":" + str(path),
                "saved_pointer": {"expression_row": row_id, "native_node": native_value, "path": path},
                "native_node": native_value,
                "path": path,
                "role": row.get("role"),
                "source": row.get("source", "physical_source_capsule"),
                "source_R_is_profile_not_physical_r": native == "241787",
            },
        )
        self.primitive_cache[row_id] = result
        return result

    def _coordinate(self, row_id: str, row: Mapping[str, Any]) -> str:
        name = str(row.get("name", ""))
        if name not in {"tau", "theta", "X_phys", "Y_phys", "z_phys", "t_phys"}:
            raise KernelInputError("unsupported free coordinate " + name)
        if name in self.coordinate_cache:
            return self.coordinate_cache[name]
        result = self.new(
            "coordinate_port",
            (),
            {
                "name": name,
                "saved_pointer": row_id,
                "role": row.get("role"),
                "physical_theta_not_native_phase": name == "theta",
            },
        )
        self.coordinate_cache[name] = result
        return result

    def compile(self, ref: Any) -> str:
        row_id = self.row_id(ref)
        if row_id in self.cache:
            return self.cache[row_id]
        row = self.saved_by_id[row_id]
        operation = row.get("operation")
        if operation not in SUPPORTED_REACHABLE_OPERATIONS:
            raise KernelInputError("unsupported reachable saved operation " + repr(operation) + " at " + row_id)
        if operation == "native_function_reference":
            result = self._primitive(row_id, row)
        elif operation == "free_variable":
            result = self._coordinate(row_id, row)
        elif operation == "bound_variable":
            name = str(row.get("name", ""))
            if name not in {"tau", "theta"}:
                raise KernelInputError("unsupported reachable bound variable " + name)
            result = self.new(
                "coordinate_port",
                (),
                {"name": name, "saved_pointer": row_id, "bound_variable": True, "physical_theta_not_native_phase": name == "theta"},
            )
        elif operation == "rational":
            numerator, denominator = row.get("numerator"), row.get("denominator", 1)
            if not isinstance(numerator, int) or not isinstance(denominator, int) or not denominator:
                raise KernelInputError("malformed saved rational at " + row_id)
            result = self.rational(numerator, denominator, premise="saved exact algebra literal")
        else:
            args = [self.compile(arg) for arg in self._args(row)]
            attrs: dict[str, Any] = {"saved_pointer": row_id}
            if operation in ANALYTIC_OPERATIONS:
                attrs["name"] = operation
            result = self.new(str(operation), args, attrs)
        self.cache[row_id] = result
        return result

    def source_port(self, operation: str, *, opaque_id: str, saved_pointer: Any, attrs: Mapping[str, Any] | None = None) -> str:
        row_attrs: dict[str, Any] = {
            "opaque_id": opaque_id,
            "saved_pointer": saved_pointer,
            "source_values_unevaluated": True,
        }
        if attrs:
            row_attrs.update(dict(attrs))
        return self.new(operation, (), row_attrs)


class _Compiler:
    def __init__(self, capsule: Mapping[str, Any], callback: Callable[[str, Mapping[str, Any]], Any] | None) -> None:
        self.capsule = _as_mapping(capsule, "capsule")
        if callback is not None and not callable(callback):
            raise KernelInputError("callback must be callable or None")
        self.callback = callback
        self.events: list[dict[str, Any]] = []
        self.internal_counts: dict[str, int] = {
            "saved_rows": 0,
            "reachable_rows": 0,
            "fresh_ir_rows": 0,
            "primitive_source_ports": 0,
            "source_partial_joins": 0,
            "reverse_equation_receipts": 0,
            "branch_gradient_entries": 0,
            "pressure_gradient_entries": 0,
            "branch_curl_entries": 0,
            "branch_integrands": 0,
            "composition_rows": 0,
            "saved_geometry_differentials_consumed": 1,
            "geometry_dual_calls": 0,
        }
        self.saved_result, self.payload = self._locate_saved_result()
        rows = self.payload.get("expression_IR")
        if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes, bytearray)) or len(rows) != 205:
            raise KernelInputError("saved payload.expression_IR must contain exactly 205 rows")
        self.builder = _IRBuilder(rows)
        self.internal_counts["saved_rows"] = len(rows)
        self.recipe_mode = "source_partial_recipes" in self.capsule
        self.recipes = self._locate_recipes()
        self.legacy_partials = self._locate_legacy_partials()
        self.fixed_ports = _as_mapping(self.capsule.get("fixed_ports", {}), "fixed_ports")
        self.frame_guard = self._frame_guard()
        self.frame_witness = self.capsule.get("source_frame_witness")
        self.geometry_row = self._locate_geometry()
        self.inverse = self._locate_inverse()
        self.branches = _as_mapping(self.payload.get("branches"), "saved branches")
        self.radial_result = self._mapping_or_empty("radial_result", "saved_radial_result")
        self.flux_payload = self._mapping_or_empty("flux_payload", "saved_flux_payload", "flux_payload_selections", "saved_flux_source_proof_rows")
        self.physical_source_capsule = self._mapping_or_empty("physical_source_capsule", "source_capsule")
        self.reverse_receipts: list[dict[str, Any]] = []
        self.branch_results: dict[str, dict[str, Any]] = {}
        self.positive_cache: dict[str, bool] = {}
        self.fraction_cache: dict[str, Fraction | None] = {}
        self.frame_certificate: dict[str, Any] | None = None
        self.frame_receipts: list[dict[str, Any]] = []

    def _locate_saved_result(self) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
        candidates = [self.capsule.get(key) for key in ("saved_exterior_result", "saved_result", "exterior_result", "saved_exterior_payload") if key in self.capsule]
        candidates.append(self.capsule)
        for candidate in candidates:
            if not _is_mapping(candidate):
                continue
            result = candidate.get("result")
            if _is_mapping(result) and _is_mapping(result.get("payload")):
                return candidate, result["payload"]
            if _is_mapping(candidate.get("payload")) and "expression_IR" in candidate["payload"]:
                return candidate, candidate["payload"]
            if "expression_IR" in candidate and "branches" in candidate:
                return candidate, candidate
        raise KernelInputError("saved exterior payload is missing")

    def _locate_recipes(self) -> Mapping[str, Any] | None:
        value = self.capsule.get("source_partial_recipes")
        if value is None:
            return None
        return _as_mapping(value, "source_partial_recipes")

    def _locate_legacy_partials(self) -> Any:
        for key in ("source_partials", "source_partial_table", "source_partials_table", "source_partial_rows", "partials"):
            if key in self.capsule:
                return self.capsule[key]
        if self.recipes is not None:
            return None
        for key in ("source_partials", "source_partial_table", "source_partials_table"):
            if key in self.payload:
                return self.payload[key]
        raise KernelInputError("source partial recipes/table is missing")

    def _mapping_or_empty(self, *keys: str) -> Mapping[str, Any]:
        for key in keys:
            if key in self.capsule:
                return _as_mapping(self.capsule[key], key)
        return {}

    def _frame_guard(self) -> Mapping[str, Any]:
        for key in ("source_frame_witness", "frame_premise", "frame_contract", "ordinary_Z_frame_premise"):
            value = self.capsule.get(key)
            if _is_mapping(value):
                return value
        return {}

    def _frame_records(self) -> Mapping[str, Any]:
        witness = _as_mapping(self.frame_witness, "source_frame_witness")
        records = witness.get("source_graph_records", witness.get("raw_graph_records"))
        if not _is_mapping(records) or not records:
            raise KernelInputError("source_frame_witness source_graph_records is missing")
        normalized: dict[str, Any] = {}
        for raw_id, row in records.items():
            try:
                key = _native_id(raw_id)
            except KernelInputError as error:
                raise KernelInputError("source-frame graph id is not numeric: " + repr(raw_id)) from error
            normalized[key] = row
        return normalized

    def _frame_root(self, roots: Mapping[str, Any], name: str, default: int) -> str:
        value = roots.get(name)
        if value is None:
            aliases = {
                "R": ("native_R", "profile_R", "R_node"),
                "actual_phase": ("phase", "actual_phase_node"),
            }.get(name, ())
            for alias in aliases:
                if alias in roots:
                    value = roots[alias]
                    break
        if value is None:
            # Accept the explicit compact labels used by handoff capsules.
            for key, candidate in roots.items():
                if _compact(key) == _compact(name) or _compact(key).startswith(_compact(name) + "node"):
                    value = candidate
                    break
        if value is None:
            value = default
        return _native_id(value)

    def _frame_atoms(self, witness: Mapping[str, Any]) -> set[str]:
        raw_atoms = witness.get("fixed_parameter_atoms", witness.get("fixed_source_atoms", ()))
        if not isinstance(raw_atoms, Sequence) or isinstance(raw_atoms, (str, bytes, bytearray)):
            raise KernelInputError("source-frame fixed_parameter_atoms is missing")
        atoms: set[str] = set()
        for item in raw_atoms:
            if _is_mapping(item):
                item = item.get("node_id", item.get("native_node", item.get("node")))
            atoms.add(_native_id(item))
        if not atoms:
            raise KernelInputError("source-frame fixed_parameter_atoms is empty")
        return atoms

    def _frame_args(self, row: Mapping[str, Any], operation: str) -> list[str]:
        if operation in {"negative", "analytic_unary"}:
            value = row.get("argument")
            if value is None:
                raise KernelInputError("source-frame " + operation + " lacks argument")
            return [_native_id(value)]
        if operation in {"sum", "product", "multiply", "add"}:
            values = row.get("arguments")
            if not isinstance(values, Sequence) or isinstance(values, (str, bytes, bytearray)):
                raise KernelInputError("source-frame " + operation + " lacks arguments")
            return [_native_id(value) for value in values]
        if operation in {"positive_quotient", "quotient"}:
            if row.get("numerator") is None or row.get("denominator") is None:
                raise KernelInputError("source-frame positive_quotient lacks numerator/denominator")
            return [_native_id(row["numerator"]), _native_id(row["denominator"])]
        if operation == "function_substitution":
            return [_native_id(row.get("expression")), _native_id(row.get("variable")), _native_id(row.get("value"))]
        return []

    def _frame_off_seam(self, row: Mapping[str, Any], witness: Mapping[str, Any]) -> bool:
        for key in ("off_phase_seam", "phase_off_seam", "regular_open_phase_cell", "phase_cell_off_seam"):
            if row.get(key) is True or witness.get(key) is True:
                return True
        phase_scope = witness.get("phase_scope")
        if _is_mapping(phase_scope) and bool(phase_scope.get("off_seam", phase_scope.get("regular_open_cell", False))):
            return True
        source_frame = self.capsule.get("source_frame")
        required_domain = source_frame.get("required_domain", "") if _is_mapping(source_frame) else ""
        return "regularopen" in _compact(required_domain)

    def _frame_zero_differential(self) -> dict[str, Any]:
        """Induct a fresh conditional D_Z=0 sequent over the frame witness.

        This is a local structural proof over the caller-supplied graph only.
        It treats fixed source atoms as conditional inputs, keeps N opaque, and
        rejects Z or any unsupported/missing graph edge rather than inferring a
        zero from metadata.
        """

        if self.frame_witness is None:
            raise KernelInputError("source_frame_witness is required before admitting spatial partials")
        witness = _as_mapping(self.frame_witness, "source_frame_witness")
        records = self._frame_records()
        roots = _as_mapping(witness.get("roots"), "source-frame roots")
        ordinary_Z = _native_id(witness.get("ordinary_Z_node", 20))
        native_x = _native_id(witness.get("native_x_node", 241536))
        if ordinary_Z != "20":
            raise KernelInputError("source-frame ordinary_Z_node must be 20")
        atoms = self._frame_atoms(witness)
        receipts: list[dict[str, Any]] = []
        cache: dict[str, str] = {}
        active: set[str] = set()
        order: list[str] = []
        substitution_nodes = witness.get("substitution_binding_nodes", ())
        if not isinstance(substitution_nodes, Sequence) or isinstance(substitution_nodes, (str, bytes, bytearray)):
            substitution_nodes = ()
        allowed_bound = {native_x, _native_id(witness.get("fraction_variable_node", 84528))}
        allowed_bound.update(_native_id(item) for item in substitution_nodes)

        def visit(node: str) -> str:
            if node == ordinary_Z:
                raise KernelInputError("source-frame zero induction encountered ordinary Z node 20")
            if node in cache:
                return cache[node]
            if node in active:
                raise KernelInputError("cycle in source-frame graph at " + node)
            if node not in records:
                raise KernelInputError("source-frame graph node is missing " + node)
            row = _as_mapping(records[node], "source-frame graph row " + node)
            operation = row.get("operation")
            if not isinstance(operation, str):
                raise KernelInputError("source-frame graph operation is missing at " + node)
            active.add(node)
            children: list[str] = []
            if operation in {"negative", "sum", "product", "multiply", "add", "positive_quotient", "quotient", "analytic_unary", "function_substitution"}:
                children = self._frame_args(row, operation)
            child_refs = [visit(child) for child in children]
            if operation == "exact_rational":
                rule = "exact rational literal has zero ordinary-Z differential"
            elif operation == "exact_positive_integer_expression":
                rule = "opaque exact positive integer expression is held fixed; N is not materialized"
            elif operation == "current_original_source_parameter":
                if node not in atoms:
                    raise KernelInputError("unlisted source parameter atom " + node)
                rule = "listed fixed source atom has conditional zero ordinary-Z differential"
            elif operation == "bound_variable":
                if node == ordinary_Z:
                    raise KernelInputError("ordinary Z bound variable is not a fixed atom")
                if node not in allowed_bound:
                    raise KernelInputError("unadmitted source-frame bound variable " + node)
                rule = "declared held native-x/substitution bound variable has zero ordinary-Z differential"
            elif operation == "negative":
                if len(child_refs) != 1:
                    raise KernelInputError("source-frame negative arity mismatch at " + node)
                rule = "negative rule preserves zero differential"
            elif operation in {"sum", "add"}:
                if not child_refs:
                    raise KernelInputError("source-frame sum is empty at " + node)
                rule = "sum rule with every child differential zero"
            elif operation in {"product", "multiply"}:
                if len(child_refs) < 2:
                    raise KernelInputError("source-frame product arity mismatch at " + node)
                rule = "product rule with every factor differential zero"
            elif operation in {"positive_quotient", "quotient"}:
                if len(child_refs) != 2:
                    raise KernelInputError("source-frame quotient arity mismatch at " + node)
                rule = "positive quotient rule with numerator/denominator differentials zero"
            elif operation == "analytic_unary":
                name = row.get("name")
                if name not in {"exp", "log", "fractional_part"}:
                    raise KernelInputError("unsupported source-frame analytic unary " + repr(name))
                if name == "fractional_part" and not self._frame_off_seam(row, witness):
                    raise KernelInputError("fractional_part frame row lacks off-seam premise at " + node)
                rule = "analytic unary " + str(name) + " chain rule with zero argument differential"
            elif operation == "function_substitution":
                if len(child_refs) != 3:
                    raise KernelInputError("source-frame function_substitution arity mismatch at " + node)
                rule = "function-substitution derivative evaluated from expression/value/zero bound-variable tangent"
            else:
                raise KernelInputError("unsupported source-frame operation " + repr(operation) + " at " + node)
            zero = self.builder.rational(0, premise="source-frame D_Z=0 induction: " + rule)
            cache[node] = zero
            active.remove(node)
            order.append(node)
            receipts.append({"node": int(node) if node.lstrip("-").isdigit() else node, "operation": operation, "children": [int(child) if child.lstrip("-").isdigit() else child for child in children], "child_zero_refs": child_refs, "zero_ref": zero, "rule": rule, "all_rule_premises_zero": True, "source_values_unevaluated": True})
            return zero

        root_ids = {
            "Pstar": self._frame_root(roots, "Pstar", 33572),
            "delta": self._frame_root(roots, "delta", 25220),
            "Rm": self._frame_root(roots, "Rm", 85743),
            "Ra": self._frame_root(roots, "Ra", 65604),
            "phase_origin": self._frame_root(roots, "phase_origin", 24954),
            "R": self._frame_root(roots, "R", 241787),
            "actual_phase": self._frame_root(roots, "actual_phase", 242551),
        }
        root_zero_refs = {name: visit(node) for name, node in root_ids.items()}
        self.frame_receipts = receipts
        self.frame_certificate = {
            "conditional": True,
            "roots": {name: int(node) if node.lstrip("-").isdigit() else node for name, node in root_ids.items()},
            "root_zero_refs": root_zero_refs,
            "ordinary_Z_node": int(ordinary_Z),
            "native_x_node": int(native_x),
            "fixed_parameter_atoms": sorted(int(atom) for atom in atoms if atom.lstrip("-").isdigit()),
            "fixed_atom_premise": witness.get("fixed_atom_premise"),
            "opaque_N_node": witness.get("opaque_N_node", 18),
            "substitution_binding_nodes": [int(item) if _native_id(item).lstrip("-").isdigit() else _native_id(item) for item in substitution_nodes],
            "reachable_node_order": [int(node) if node.lstrip("-").isdigit() else node for node in order],
            "receipts": receipts,
            "conditional_D_Z_zero": {name: True for name in root_ids},
            "conditional_Rm_Ra_phase_R_invariance": True,
            "conditional_fixed_x_Z_equals_fixed_R_Z": True,
            "connective_equation": "R=Rm*exp(17+exp(40)+native_x); D_x R=R>0 on the source domain. At fixed R, native_x=log(R/Rm)-17-exp(40) gives D_Z native_x=-D_Z Rm/Rm=0. Hence ordinary D_Z at fixed x agrees with D_Z at fixed R under the proved local fixed-parameter frame.",
            "fractional_part_off_seam_only": True,
            "N_materialized": False,
            "unconditional_Rm_Z_certificate": False,
            "unconditional_source_Z_admission": False,
            "global_phase_seam_regularity": False,
            "source_values_unevaluated": True,
        }
        return self.frame_certificate

    def _locate_geometry(self) -> Mapping[str, Any]:
        direct = self.capsule.get("geometry_dual_row")
        if _is_mapping(direct):
            return direct
        callbacks = self.saved_result.get("callback_events") if _is_mapping(self.saved_result) else None
        if isinstance(callbacks, Sequence) and not isinstance(callbacks, (str, bytes, bytearray)):
            for event in callbacks:
                if _is_mapping(event) and event.get("event") == "geometry_dual" and _is_mapping(event.get("row")):
                    return event["row"]
        direct = self.payload.get("geometry_dual")
        if _is_mapping(direct):
            return direct
        raise KernelInputError("geometry_dual row is missing")

    def _locate_inverse(self) -> Mapping[str, Any]:
        inverse = _as_mapping(self.payload.get("inverse"), "saved inverse")
        if inverse.get("bound_root") is None or not _is_mapping(inverse.get("bindings")):
            raise KernelInputError("saved inverse lacks explicit bindings and bound_root")
        return inverse

    def _emit(self, event: str, row: Mapping[str, Any]) -> None:
        if event not in CALLBACK_EVENT_COUNTS:
            raise KernelInputError("event outside frozen callback contract " + event)
        item = {"event": event, "row": dict(row)}
        self.events.append(item)
        if self.callback is not None:
            self.callback(event, item["row"])

    def _event_counts(self) -> dict[str, int]:
        return {event: sum(1 for item in self.events if item.get("event") == event) for event in CALLBACK_EVENT_COUNTS}

    def _manifest_counts(self, observed: Mapping[str, int]) -> dict[str, int]:
        derived = {name: 0 for name in SUCCESS_COUNTS}
        for event, targets in COUNT_EVENT_MAP.items():
            for target in targets:
                derived[target] += int(observed.get(event, 0))
        for name in SUCCESS_COUNTS:
            if SUCCESS_COUNTS[name] == 0:
                derived[name] = 0
        if derived != SUCCESS_COUNTS:
            raise KernelInputError("callback counters do not match frozen contract")
        return dict(SUCCESS_COUNTS)

    def _partial_entry_legacy(self, node: str, path: Any, branch: str) -> Mapping[str, Any] | None:
        table = self.legacy_partials

        def direct(container: Any, key: Any) -> Any:
            if _is_mapping(container):
                if key in container:
                    return container[key]
                text = str(key)
                if text in container:
                    return container[text]
            return None

        def scan(container: Any) -> Any:
            if not isinstance(container, Sequence) or isinstance(container, (str, bytes, bytearray)):
                return None
            for item in container:
                if not _is_mapping(item):
                    continue
                item_node = item.get("native_node", item.get("node", item.get("source_node")))
                if item_node is not None and _native_id(item_node) == node and item.get("branch") in (None, branch):
                    return item
            return None

        containers = [table]
        if _is_mapping(table):
            for key in ("by_native_node", "native_nodes", "ports", "fields", "entries", "rows", branch):
                if key in table:
                    containers.append(table[key])
        for container in containers:
            hit = direct(container, node) or scan(container)
            if hit is not None:
                return _as_mapping(hit, "legacy source partial entry")
            if _is_mapping(container):
                path_hit = direct(container, path)
                if path_hit is not None:
                    return _as_mapping(path_hit, "legacy source partial path entry")
        return None

    def _axis_value(self, entry: Mapping[str, Any], axis: str) -> Any:
        for nested_key in ("partials", "derivatives", "source_partials", "values"):
            nested = entry.get(nested_key)
            if _is_mapping(nested):
                for key in ("d" + axis, axis, "partial_" + axis, "profile_" + axis, "D_" + axis):
                    if key in nested:
                        return nested[key]
        for key in ("d" + axis, axis, "partial_" + axis, "profile_" + axis, "D_" + axis):
            if key in entry:
                return entry[key]
        return None

    def _frame_ok(self, recipe: Mapping[str, Any]) -> bool:
        source_frame = recipe.get("source_frame")
        if source_frame != FRAME_NAME:
            return False
        # The recipe's named frame is a declared conditional premise.  A
        # separate witness may strengthen it, but is not silently inferred.
        return bool(recipe.get("frame_premise", True))

    def _fixed_parameter_frame_ok(self, recipe: Mapping[str, Any]) -> bool:
        """Require the explicit conditional fixed-parameter frame contract."""

        admission = recipe.get("admission", recipe.get("frame_admission"))
        if admission not in {
            "declared_fixed_parameter_frame_only",
            "conditional_on_new_same_source_zero_differential_certificate",
        }:
            return False
        source_frame = self.capsule.get("source_frame")
        if not _is_mapping(source_frame):
            return False
        return source_frame.get("admission") in {
            "conditional_declared_fixed_parameter_frame",
            "conditional_declared_source_jet_frame_and_new_same_source_zero_differential",
        }

    def _pointer_lookup(self, pointer: Any) -> Any:
        """Resolve a declarative saved-flux pointer inside literal data."""

        if not isinstance(pointer, str):
            return None
        pieces = [piece for piece in pointer.split(".") if piece]
        while pieces and pieces[0] in {"result", "payload", "source_proof_rows"}:
            pieces.pop(0)
        current: Any = self.flux_payload
        for piece in pieces:
            if isinstance(current, Mapping) and piece in current:
                current = current[piece]
            elif isinstance(current, Sequence) and not isinstance(current, (str, bytes, bytearray)) and piece.isdigit() and int(piece) < len(current):
                current = current[int(piece)]
            else:
                return None
        return current

    def _partial_literal(self, value: Any, *, node: str, axis: str, role: str) -> str:
        fraction = _exact_fraction(value)
        if fraction is not None:
            return self.builder.rational(fraction.numerator, fraction.denominator, premise=role)
        if isinstance(value, Mapping):
            for key in ("saved_expression_ref", "expression_ref", "saved_pointer", "ref"):
                if key in value and _saved_id(value[key]) in self.builder.saved_by_id:
                    return self.builder.compile(value[key])
            if value.get("kind") in {"fixed", "literal"} and "value" in value:
                return self._partial_literal(value["value"], node=node, axis=axis, role=role)
            return self.builder.source_port(
                "source_partial_port",
                opaque_id="source_partial:" + node + ":" + axis + ":" + str(value),
                saved_pointer=value,
                attrs={"native_node": int(node) if node.lstrip("-").isdigit() else node, "axis": axis, "lineage": role},
            )
        saved_ref = _saved_id(value)
        if saved_ref in self.builder.saved_by_id:
            return self.builder.compile(saved_ref)
        return self.builder.source_port(
            "source_partial_port",
            opaque_id="source_partial:" + node + ":" + axis + ":" + str(value),
            saved_pointer=value,
            attrs={"native_node": int(node) if node.lstrip("-").isdigit() else node, "axis": axis, "lineage": role},
        )

    def _source_recipe(self, node: str) -> Mapping[str, Any]:
        if self.recipes is None:
            raise KernelInputError("source recipes are unavailable")
        value = self.recipes.get(node, self.recipes.get(int(node)) if node.lstrip("-").isdigit() else None)
        if value is None:
            raise KernelInputError("unknown source partial recipe for native node " + node)
        recipe = _as_mapping(value, "source partial recipe")
        if recipe.get("native_node") is not None and _native_id(recipe.get("native_node")) != node:
            raise KernelInputError("source partial recipe native_node disagrees with its key " + node)
        return recipe

    def _radial_target(self, target_name: str) -> Mapping[str, Any]:
        payload = self.radial_result.get("result", {}).get("payload", {}) if _is_mapping(self.radial_result.get("result")) else self.radial_result.get("payload", self.radial_result)
        if not _is_mapping(payload):
            raise KernelInputError("saved radial result payload is missing")
        targets = payload.get("targets", [])
        if not isinstance(targets, Sequence):
            raise KernelInputError("saved radial result targets are missing")
        for target in targets:
            if _is_mapping(target) and target.get("name") == target_name:
                tangent = target.get("tangent_ref")
                if tangent is None:
                    raise KernelInputError("saved radial target lacks tangent_ref " + target_name)
                return {
                    "target_name": target_name,
                    "tangent_ref": tangent,
                    "saved_pointer": "radial_result.result.payload.targets." + target_name,
                    "full_ir": payload.get("tangent_ir", []),
                }
        raise KernelInputError("saved radial target not found " + target_name)

    def _flux_expression(self, recipe: Mapping[str, Any]) -> Mapping[str, Any]:
        """Resolve one exact saved source row; no first-match fallback."""
        branch = recipe.get("branch")
        field = recipe.get("field")
        order = recipe.get("ordinary_Z_order")
        if branch not in {"background", "modified"} or not isinstance(field, str) or not isinstance(order, int):
            raise KernelInputError("saved flux recipe lacks exact branch/field/order")
        pointer = recipe.get("saved_result_pointer", recipe.get("saved_pointer"))
        expected = "result.payload.source_proof_rows." + branch + "." + field + "." + str(order) + ".source_expression"
        if pointer != expected:
            raise KernelInputError("saved flux pointer does not match its exact branch/field/order")
        pointed = self._pointer_lookup(pointer)
        if not _is_mapping(pointed):
            raise KernelInputError("exact saved flux source expression is missing: " + expected)
        return {"saved_pointer": pointer, "source_expression": pointed, "branch": branch, "field": field, "ordinary_Z_order": order, "selection": "exact branch/field/ordinary_Z_order; no fallback"}

    def _resolve_source_field(self, recipe: Mapping[str, Any], axis: str, branch: str, node: str) -> str:
        if not self._frame_ok(recipe):
            raise KernelInputError("source field lacks declared fixed-profile-R frame premise")
        if axis == "theta":
            theta = recipe.get("dtheta")
            if theta != "typed_zero" and not (_is_mapping(theta) and theta.get("kind") == "typed_zero"):
                raise KernelInputError("source theta partial is not an explicit typed_zero")
            return self.builder.rational(0, premise="typed_zero source theta partial; source value retained opaque")
        if axis == "R":
            radial = _as_mapping(recipe.get("radial"), "source field radial recipe")
            kind = radial.get("kind")
            if kind == "saved_radial_variation":
                target_name = radial.get("target_name")
                if not isinstance(target_name, str):
                    raise KernelInputError("saved_radial_variation target_name is missing")
                target = self._radial_target(target_name)
                target["saved_result_sha256"] = radial.get("saved_result_sha256")
            elif kind == "saved_flux_source_expression":
                target = self._flux_expression({**recipe, **radial})
            else:
                raise KernelInputError("unknown source radial recipe kind " + repr(kind))
            dy = self.builder.source_port(
                "source_radial_variation_port",
                opaque_id="radial_variation:" + node + ":" + str(target.get("saved_pointer")),
                saved_pointer=target.get("saved_pointer"),
                attrs={"target_name": target.get("target_name"), "tangent_ref": target.get("tangent_ref"), "full_ir_data": True, "saved_result_sha256": target.get("saved_result_sha256"), "source_expression_data": target.get("source_expression"), "source_frame": FRAME_NAME, "radial_equation": ("E_R=newDyE/R" if recipe.get("field") == "E" else "V_R=newDyV/R" if recipe.get("field") == "V" else "m_R=(V-m)/R" if recipe.get("field") == "m" and recipe.get("ordinary_Z_order") == 0 else "m_Z1_R=(V_Z1-m_Z1)/R if commutation admitted" if recipe.get("field") == "m" and recipe.get("ordinary_Z_order") == 1 else "p_R=accepted_FTOC_source_expression/R"), "commutation_premise": recipe.get("ordinary_Z_order") == 1 and recipe.get("field") == "m"},
            )
            R_value = self._profile_R_value()
            return self.builder.new("quotient", [dy, R_value], {"equation": "D_R=D_y/R", "radial_variation_data": True, "source_frame": FRAME_NAME})
        axial = _as_mapping(recipe.get("axial"), "source field axial recipe")
        if axial.get("kind") != "native_ordinary_Z_successor":
            raise KernelInputError("unknown source axial recipe kind " + repr(axial.get("kind")))
        successor = axial.get("node")
        if successor is None:
            raise KernelInputError("ordinary Z successor node is missing")
        return self.builder.source_port(
            "source_axial_successor_port",
            opaque_id="ordinary_Z_successor:" + node + ":" + str(successor),
            saved_pointer={"native_node": successor, "branch": branch, "field": recipe.get("field"), "ordinary_Z_order": recipe.get("ordinary_Z_order")},
            attrs={"native_node": int(_native_id(successor)) if _native_id(successor).lstrip("-").isdigit() else successor, "source_frame": FRAME_NAME, "ordinary_Z_successor": True, "source_frame_admitted_unconditionally": False},
        )

    def _resolve_alias(self, recipe: Mapping[str, Any], axis: str, branch: str, node: str) -> str:
        if recipe.get("kind") != "accepted_my_alias":
            raise KernelInputError("unknown alias recipe kind")
        source_frame = self.capsule.get("source_frame")
        alias_frame_ok = _is_mapping(source_frame) and source_frame.get("admission") in {
            "conditional_declared_fixed_parameter_frame",
            "conditional_declared_source_jet_frame_and_new_same_source_zero_differential",
        }
        if not alias_frame_ok:
            raise KernelInputError("accepted_my_alias lacks declared frame premise")
        V_node = str(recipe.get("V_node", 244378))
        m_node = str(recipe.get("m_node", 244710))
        V_recipe = self._source_recipe(V_node)
        m_recipe = self._source_recipe(m_node)
        V_ref = self._resolve_recipe(V_recipe, axis, branch, V_node)
        m_ref = self._resolve_recipe(m_recipe, axis, branch, m_node)
        return self.builder.new("add", [V_ref, self.builder.new("negative", [m_ref], {"equation": "-modified_m0_partial"})], {"equation": "accepted_my=modified_V0-modified_m0", "accepted_alias": True, "new_subtraction_inside_kernel": True, "source_frame": FRAME_NAME})

    def _resolve_recipe(self, recipe: Mapping[str, Any], axis: str, branch: str, node: str) -> str:
        kind = recipe.get("kind")
        if kind == "coordinate":
            named_axis = recipe.get("axis")
            value = recipe.get("d" + axis, recipe.get(axis))
            if value is None:
                value = 1 if named_axis == axis else 0
            return self._partial_literal(value, node=node, axis=axis, role="explicit coordinate recipe")
        if kind == "fixed_parameter":
            if not self._fixed_parameter_frame_ok(recipe):
                raise KernelInputError("fixed parameter recipe lacks explicit frame guard")
            value = recipe.get("d" + axis, recipe.get(axis))
            if value is None:
                value = 0
            return self._partial_literal(value, node=node, axis=axis, role="explicit fixed parameter recipe")
        if kind == "source_exact_constant":
            raw = _as_mapping(recipe.get("raw_record"), "source exact constant raw_record")
            fraction = _exact_fraction(raw)
            if fraction is None:
                raise KernelInputError("source exact constant lacks exact rational raw_record")
            return self._partial_literal(0, node=node, axis=axis, role="explicit source exact constant partial")
        if kind == "log_profile_coordinate":
            raise KernelInputError("native_x derivative is outside the admitted physical branch partial contract")
        if kind == "source_field":
            return self._resolve_source_field(recipe, axis, branch, node)
        if kind == "accepted_my_alias":
            return self._resolve_alias(recipe, axis, branch, node)
        raise KernelInputError("unknown source partial recipe kind " + repr(kind))

    def _native_partial(self, row: Mapping[str, Any], axis: str, branch: str) -> str:
        node = _native_id(row.get("native_node"))
        if self.recipes is not None:
            recipe = self._source_recipe(node)
            result = self._resolve_recipe(recipe, axis, branch, node)
        else:
            path = row.get("path")
            entry = self._partial_entry_legacy(node, path, branch)
            value = self._axis_value(entry, axis) if entry is not None else None
            if value is None and node in FIXED_PORT_PARTIALS:
                value = FIXED_PORT_PARTIALS[node][axis]
            if value is None:
                override = self.fixed_ports.get(node, self.fixed_ports.get(int(node)) if node.lstrip("-").isdigit() else None)
                if _is_mapping(override):
                    value = override.get("d" + axis, override.get(axis))
            if value is None:
                raise KernelInputError("unknown source partial for native node " + node + " axis " + axis)
            result = self._partial_literal(value, node=node, axis=axis, role="legacy supplied source partial")
        self.internal_counts["source_partial_joins"] += 1
        return result

    def _coordinate_partial(self, row: Mapping[str, Any], axis: str) -> str:
        name = str(row.get("name", ""))
        if name == "theta":
            value = 1 if axis == "theta" else 0
        elif name in {"tau", "X_phys", "Y_phys", "z_phys", "t_phys"}:
            value = 0
        else:
            raise KernelInputError("unknown coordinate partial " + name)
        return self.builder.rational(value, premise="fixed coordinate partial")

    def _profile_R_value(self) -> str:
        for row_id, row in self.builder.saved_by_id.items():
            if row.get("operation") == "native_function_reference" and _native_id(row.get("native_node")) == "241787":
                return self.builder.compile(row_id)
        raise KernelInputError("saved profile-R primitive is missing")

    def _theta_value(self) -> str:
        for row_id, row in self.builder.saved_by_id.items():
            if row.get("operation") == "free_variable" and row.get("name") == "theta":
                return self.builder.compile(row_id)
        raise KernelInputError("saved physical theta free variable is missing")

    def _saved_fraction(self, row_id: str) -> Fraction | None:
        if row_id in self.fraction_cache:
            return self.fraction_cache[row_id]
        row = self.builder.saved_by_id[row_id]
        operation = row.get("operation")
        result: Fraction | None
        if operation == "rational":
            result = _exact_fraction(row)
        elif operation == "negative":
            args = self.builder._args(row)
            value = self._saved_fraction(args[0]) if len(args) == 1 else None
            result = -value if value is not None else None
        elif operation in {"add", "multiply"}:
            values = [self._saved_fraction(arg) for arg in self.builder._args(row)]
            if any(value is None for value in values):
                result = None
            elif operation == "add":
                result = sum(values, Fraction(0, 1))
            else:
                result = Fraction(1, 1)
                for value in values:
                    result *= value
        else:
            result = None
        self.fraction_cache[row_id] = result
        return result

    def _is_positive(self, row_id: str) -> bool:
        if row_id in self.positive_cache:
            return self.positive_cache[row_id]
        # The saved geometry callback carries the accepted positive-Jacobian
        # domain premise even though the raw algebraic row is not syntactically
        # a sum of positive factors.
        saved_J2 = _saved_id(self.geometry_row.get("J2"))
        if saved_J2 is not None and row_id == saved_J2:
            self.positive_cache[row_id] = True
            return True
        row = self.builder.saved_by_id[row_id]
        operation = row.get("operation")
        result = False
        if operation == "rational":
            value = _exact_fraction(row)
            result = value is not None and value > 0
        elif operation == "native_function_reference":
            node = _native_id(row.get("native_node"))
            path = _compact(row.get("path"))
            result = node in {"241787", "33572", "85743"} or any(tag in path for tag in ("positive_r", "pstar"))
        elif operation == "free_variable":
            result = row.get("name") == "tau"
        elif operation == "bound_variable":
            result = row.get("name") == "tau"
        elif operation in {"multiply", "quotient"}:
            args = self.builder._args(row)
            result = bool(args) and all(self._is_positive(arg) for arg in args)
        elif operation == "power":
            args = self.builder._args(row)
            result = len(args) == 2 and self._is_positive(args[0])
        elif operation == "add":
            args = self.builder._args(row)
            if args and all(self._is_positive(arg) for arg in args):
                result = True
            elif len(args) == 2 and self._saved_fraction(args[0]) == 1:
                neg = self.builder.saved_by_id[args[1]]
                if neg.get("operation") == "negative":
                    body = self.builder._args(neg)
                    if len(body) == 1:
                        factors = self.builder._args(self.builder.saved_by_id[body[0]])
                        z_count = sum(1 for item in factors if self._is_z_ref(item))
                        result = z_count == 2 and (len(factors) == 2 or (len(factors) == 3 and any(self._is_delta_ref(item) for item in factors)))
        elif operation == "implicit_positive_root":
            result = True
        self.positive_cache[row_id] = bool(result)
        return bool(result)

    def _is_z_ref(self, row_id: str) -> bool:
        row = self.builder.saved_by_id[row_id]
        return row.get("operation") == "native_function_reference" and _native_id(row.get("native_node")) == "20"

    def _is_delta_ref(self, row_id: str) -> bool:
        row = self.builder.saved_by_id[row_id]
        return row.get("operation") == "native_function_reference" and _native_id(row.get("native_node")) == "25220"

    def _reachable(self, root_id: str) -> list[str]:
        active: set[str] = set()
        seen: set[str] = set()
        post: list[str] = []

        def visit(row_id: str) -> None:
            if row_id in active:
                raise KernelInputError("cycle in saved value DAG at " + row_id)
            if row_id in seen:
                return
            row = self.builder.saved_by_id.get(row_id)
            if row is None:
                raise KernelInputError("saved DAG references unknown row " + row_id)
            operation = row.get("operation")
            if operation not in SUPPORTED_REACHABLE_OPERATIONS:
                raise KernelInputError("unsupported reachable operation " + repr(operation) + " at " + row_id)
            active.add(row_id)
            current = self.builder.saved_order[row_id]
            for arg in self.builder._args(row):
                if arg not in self.builder.saved_by_id:
                    raise KernelInputError("saved DAG argument is unknown " + arg)
                if self.builder.saved_order[arg] >= current:
                    raise KernelInputError("saved expression IR is not topologically ordered at " + row_id)
                visit(arg)
            active.remove(row_id)
            seen.add(row_id)
            post.append(row_id)

        visit(root_id)
        self.internal_counts["reachable_rows"] += len(post)
        return post

    def _sum_weight(self, weights: dict[str, str], row_id: str, term: str) -> None:
        if row_id in weights:
            weights[row_id] = self.builder.new("add", [weights[row_id], term], {"reverse_weight_sum": True})
        else:
            weights[row_id] = term

    def _product(self, refs: Sequence[str]) -> str:
        if not refs:
            return self.builder.rational(1, premise="empty reverse product")
        if len(refs) == 1:
            return refs[0]
        return self.builder.new("multiply", refs, {"reverse_product": True})

    def _reverse(self, root_ref: Any, branch: str, field: str) -> dict[str, Any]:
        root_id = self.builder.row_id(root_ref)
        root_row = self.builder.saved_by_id[root_id]
        if root_row.get("operation") in {"vector", "compose_native_function_and_geometry"}:
            raise KernelInputError("reverse target must be an uncomposed scalar root")
        reachable = self._reachable(root_id)
        values = {row_id: self.builder.compile(row_id) for row_id in reachable}
        weights: dict[str, str] = {root_id: self.builder.rational(1, premise="reverse target weight")}
        receipts: list[dict[str, Any]] = []

        for row_id in reversed(reachable):
            row = self.builder.saved_by_id[row_id]
            weight = weights.get(row_id)
            if weight is None:
                continue
            operation = row.get("operation")
            args = self.builder._args(row)
            propagated: dict[str, str] = {}
            propagated_contributions: dict[str, list[str]] = {}

            def note_propagation(argument: str, term: str) -> None:
                """Retain every operand contribution, including duplicate operands."""
                propagated_contributions.setdefault(argument, []).append(term)
                # Keep the historical single-value view for consumers that only
                # need one representative; the complete evidence lives beside it.
                propagated[argument] = term

            if operation in {"rational", "native_function_reference", "free_variable", "bound_variable"}:
                pass
            elif operation == "add":
                for arg in args:
                    self._sum_weight(weights, arg, weight)
                    note_propagation(arg, weight)
            elif operation == "negative":
                if len(args) != 1:
                    raise KernelInputError("negative row must have one argument")
                term = self.builder.new("negative", [weight], {"reverse_rule": "negative"})
                self._sum_weight(weights, args[0], term)
                note_propagation(args[0], term)
            elif operation == "multiply":
                if len(args) < 2:
                    raise KernelInputError("multiply row has fewer than two arguments")
                for index, arg in enumerate(args):
                    factor = self._product([values[item] for j, item in enumerate(args) if j != index])
                    term = self.builder.new("multiply", [weight, factor], {"reverse_rule": "product"})
                    self._sum_weight(weights, arg, term)
                    note_propagation(arg, term)
            elif operation == "quotient":
                if len(args) != 2:
                    raise KernelInputError("quotient row must have two arguments")
                numerator, denominator = args
                if not self._is_positive(denominator):
                    raise KernelInputError("quotient denominator lacks positive domain premise at " + row_id)
                denominator_value, numerator_value = values[denominator], values[numerator]
                num_term = self.builder.new("quotient", [weight, denominator_value], {"reverse_rule": "quotient_numerator"})
                square = self.builder.new("multiply", [denominator_value, denominator_value], {"reverse_rule": "denominator_square"})
                den_num = self.builder.new("multiply", [weight, numerator_value], {"reverse_rule": "quotient_denominator"})
                den_term = self.builder.new("negative", [self.builder.new("quotient", [den_num, square], {"reverse_rule": "quotient_denominator_fraction"})], {"reverse_rule": "quotient_denominator_sign"})
                self._sum_weight(weights, numerator, num_term)
                self._sum_weight(weights, denominator, den_term)
                note_propagation(numerator, num_term)
                note_propagation(denominator, den_term)
            elif operation == "power":
                if len(args) != 2:
                    raise KernelInputError("power row must have base and exponent")
                base, exponent = args
                exponent_value, base_value = values[exponent], values[base]
                exponent_fraction = self._saved_fraction(exponent)
                if exponent_fraction is not None:
                    if exponent_fraction != 0:
                        if (exponent_fraction.denominator != 1 or exponent_fraction < 0) and not self._is_positive(base):
                            raise KernelInputError("power base lacks positive premise at " + row_id)
                        minus_one = exponent_fraction - 1
                        power = self.builder.new("power", [base_value, self.builder.rational(minus_one.numerator, minus_one.denominator, premise="exact algebra power-chain exponent")], {"reverse_rule": "exact_rational_power_chain", "constant_exponent": str(exponent_fraction), "base_positive_or_nonzero": self._is_positive(base)})
                        factor = self.builder.new("multiply", [exponent_value, power], {"reverse_rule": "power_chain_factor"})
                        term = self.builder.new("multiply", [weight, factor], {"reverse_rule": "power_chain"})
                        self._sum_weight(weights, base, term)
                        note_propagation(base, term)
                else:
                    if not self._is_positive(base):
                        raise KernelInputError("general power base lacks positive premise at " + row_id)
                    value = values[row_id]
                    log_base = self.builder.new("log", [base_value], {"reverse_rule": "general_power_log", "positive_base": True})
                    base_factor = self.builder.new("quotient", [self.builder.new("multiply", [value, exponent_value], {"reverse_rule": "general_power_base_num"}), base_value], {"reverse_rule": "general_power_base_factor"})
                    exponent_factor = self.builder.new("multiply", [value, log_base], {"reverse_rule": "general_power_exponent_factor"})
                    base_term = self.builder.new("multiply", [weight, base_factor], {"reverse_rule": "general_power_base"})
                    exponent_term = self.builder.new("multiply", [weight, exponent_factor], {"reverse_rule": "general_power_exponent"})
                    self._sum_weight(weights, base, base_term)
                    self._sum_weight(weights, exponent, exponent_term)
                    note_propagation(base, base_term)
                    note_propagation(exponent, exponent_term)
            elif operation in ANALYTIC_OPERATIONS:
                if len(args) != 1:
                    raise KernelInputError(operation + " row must have one argument")
                argument = args[0]
                argument_value = values[argument]
                if operation == "exp":
                    factor = values[row_id]
                elif operation == "log":
                    if not self._is_positive(argument):
                        raise KernelInputError("log argument lacks positive premise at " + row_id)
                    factor = self.builder.new("quotient", [self.builder.rational(1), argument_value], {"reverse_rule": "log"})
                elif operation == "cos":
                    factor = self.builder.new("negative", [self.builder.new("sin", [argument_value], {"reverse_rule": "cos_companion"})], {"reverse_rule": "cos_sign"})
                else:
                    factor = self.builder.new("cos", [argument_value], {"reverse_rule": "sin_companion"})
                term = self.builder.new("multiply", [weight, factor], {"reverse_rule": operation})
                self._sum_weight(weights, argument, term)
                note_propagation(argument, term)
            else:
                raise KernelInputError("unsupported reverse operation " + repr(operation))
            receipts.append({"branch": branch, "field": field, "saved_node": row_id, "operation": operation, "incoming_weight": weight, "propagated_weights": propagated, "propagated_weight_contributions": {key: list(values) for key, values in propagated_contributions.items()}, "aggregated_weights_after_step": {key: weights[key] for key in propagated_contributions if key in weights}, "topological_position": self.builder.saved_order[row_id], "proof": "reverse adjoint rule over saved value DAG", "source_values_unevaluated": True})

        primitive_terms: dict[str, dict[str, str]] = {}
        for row_id in reachable:
            weight = weights.get(row_id)
            if weight is None:
                continue
            row = self.builder.saved_by_id[row_id]
            operation = row.get("operation")
            if operation == "native_function_reference":
                primitive_terms[row_id] = {axis: self.builder.new("multiply", [weight, self._native_partial(row, axis, branch)], {"reverse_rule": "primitive_port_chain", "saved_primitive": row_id, "axis": axis}) for axis in AXES}
            elif operation in {"free_variable", "bound_variable"}:
                primitive_terms[row_id] = {axis: self.builder.new("multiply", [weight, self._coordinate_partial(row, axis)], {"reverse_rule": "coordinate_port_chain", "saved_primitive": row_id, "axis": axis}) for axis in AXES}
        derivatives: dict[str, str] = {}
        for axis in AXES:
            terms = [value[axis] for value in primitive_terms.values() if axis in value]
            if not terms:
                raise KernelInputError("reverse adjoint produced no terms for " + axis)
            derivatives[axis] = terms[0] if len(terms) == 1 else self.builder.new("add", terms, {"reverse_rule": "primitive_port_sum", "axis": axis})
        self.reverse_receipts.extend(receipts)
        self.internal_counts["reverse_equation_receipts"] += len(receipts)
        native_groups: dict[str, list[str]] = {}
        for saved_primitive in primitive_terms:
            native = _native_id(self.builder.saved_by_id[saved_primitive].get("native_node"))
            native_groups.setdefault(native, []).append(saved_primitive)
        return {"field": field, "saved_value_ref": root_id, "value_ref": values[root_id], "profile_partials": derivatives, "primitive_contributions": primitive_terms, "primitive_contributions_by_native_node": native_groups, "shared_native_ids_accumulated": True, "reachable_saved_nodes": reachable, "equation_receipts": receipts, "source_values_unevaluated": True}

    def _geometry(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for name in ("r", "z", "lambda"):
            row = _as_mapping(self.geometry_row.get(name), "geometry " + name)
            result[name] = {"value": self.builder.compile(row.get("value")), "saved_value": self.builder.row_id(row.get("value")), "D_R": self.builder.compile(row.get("D_R")), "D_Z": self.builder.compile(row.get("D_Z")), "saved_D_R": self.builder.row_id(row.get("D_R")), "saved_D_Z": self.builder.row_id(row.get("D_Z"))}
        J2 = self.geometry_row.get("J2")
        result["J2"] = {"value": self.builder.compile(J2), "saved_value": self.builder.row_id(J2)}
        result["inverse"] = {"saved_pointer": "result.payload.inverse", "bound_root": self.inverse.get("bound_root"), "bound_variable": self.inverse.get("bound_variable"), "bindings": dict(self.inverse.get("bindings", {})), "native_expression_substitution": bool(self.inverse.get("native_expression_substitution", False)), "numerical_inverse_available": bool(self.inverse.get("numerical_inverse_available", False)), "positive_root_data": True, "inverse_not_differentiated": True}
        return result

    def _transform(self, profile: Mapping[str, str], geometry: Mapping[str, Any], theta: str) -> dict[str, str]:
        F_R, F_Z, F_theta = profile["R"], profile["Z"], profile["theta"]
        r, J2 = geometry["r"]["value"], geometry["J2"]["value"]
        r_R, r_Z = geometry["r"]["D_R"], geometry["r"]["D_Z"]
        z_R, z_Z = geometry["z"]["D_R"], geometry["z"]["D_Z"]
        radial = self.builder.new("quotient", [self.builder.new("add", [self.builder.new("multiply", [z_Z, F_R]), self.builder.new("negative", [self.builder.new("multiply", [z_R, F_Z])])], {"equation": "z_Z*F_R-z_R*F_Z"}), J2], {"equation": "F_r=(z_Z*F_R-z_R*F_Z)/J2", "saved_geometry_differentials": True})
        axial = self.builder.new("quotient", [self.builder.new("add", [self.builder.new("multiply", [r_R, F_Z]), self.builder.new("negative", [self.builder.new("multiply", [r_Z, F_R])])], {"equation": "r_R*F_Z-r_Z*F_R"}), J2], {"equation": "F_z=(-r_Z*F_R+r_R*F_Z)/J2", "saved_geometry_differentials": True})
        c, s = self.builder.new("cos", [theta], {"physical_theta": True}), self.builder.new("sin", [theta], {"physical_theta": True})
        s_over_r, c_over_r = self.builder.new("quotient", [s, r], {"equation": "sin(theta)/r", "r_positive_domain": True}), self.builder.new("quotient", [c, r], {"equation": "cos(theta)/r", "r_positive_domain": True})
        X = self.builder.new("add", [self.builder.new("multiply", [c, radial]), self.builder.new("negative", [self.builder.new("multiply", [s_over_r, F_theta])])], {"equation": "F_X=cos(theta)*F_r-sin(theta)*F_theta/r", "physical_theta_not_native_phase": True})
        Y = self.builder.new("add", [self.builder.new("multiply", [s, radial]), self.builder.new("multiply", [c_over_r, F_theta])], {"equation": "F_Y=sin(theta)*F_r+cos(theta)*F_theta/r", "physical_theta_not_native_phase": True})
        return {"X": X, "Y": Y, "z": axial, "r": radial}

    def _composition_binding_contract(self) -> dict[str, Any]:
        bindings = dict(self.inverse.get("bindings", {}))
        native_bindings = {
            "241787": {"saved_expression_ref": "physical_expr:3", "inverse_ref": bindings.get("physical_expr:3"), "role": "profile_R"},
            "20": {"saved_expression_ref": "physical_expr:5", "inverse_ref": bindings.get("physical_expr:5"), "role": "ordinary_Z"},
            "241536": {"saved_expression_ref": "physical_expr:7", "inverse_ref": bindings.get("physical_expr:7"), "role": "native_x"},
            "physical_tau": {"saved_expression_ref": "physical_expr:8", "inverse_ref": bindings.get("physical_expr:8"), "role": "physical_time_remaining"},
            "physical_theta": {"saved_expression_ref": "physical_expr:9", "inverse_ref": bindings.get("physical_expr:9"), "role": "physical_azimuth_not_native_phase"},
        }
        source_rows = []
        for row in self.builder.rows:
            if row.get("operation") in {"primitive_source_port", "source_radial_variation_port", "source_axial_successor_port", "source_partial_port"}:
                source_rows.append({"fresh_row": row["id"], "saved_pointer": row.get("saved_pointer"), "native_node": row.get("native_node"), "opaque_id": row.get("opaque_id"), "recursive_binding": True})
        return {
            "native_coordinate_bindings": native_bindings,
            "explicit_native_bindings": native_bindings,
            "saved_inverse_bindings": bindings,
            "inverse_binding_application": "physical_expr:3/5/7/8/9 -> saved inverse substitutions; apply recursively inside source-function definitions",
            "source_function_definitions_rebound": source_rows,
            "recursive_source_function_binding": True,
            "native_references_are_opaque_functions_not_numeric_constants": True,
            "substitution_applies_inside_radial_and_axial_source_definitions": True,
        }

    def _compose(self, label: str, scalar: str, saved_ref: Any, inverse_args: Sequence[Any]) -> str:
        arguments = [self.builder.compile(item) for item in inverse_args]
        binding_contract = self._composition_binding_contract()
        return self.builder.new("compose_native_function_and_geometry", [scalar], {"label": label, "source_expression_ref": _json_scalar(saved_ref), "physical_arguments": arguments, "saved_inverse_pointer": "result.payload.inverse", "saved_positive_inverse": {"bound_root": self.inverse.get("bound_root"), "bound_variable": self.inverse.get("bound_variable"), "bindings": dict(self.inverse.get("bindings", {})), "native_expression_substitution": bool(self.inverse.get("native_expression_substitution", False)), "numerical_inverse_available": bool(self.inverse.get("numerical_inverse_available", False))}, "explicit_bindings": True, "binding_contract": binding_contract, "inverse_is_data_only": True, "circular": False})

    def _branch(self, branch: str, geometry: Mapping[str, Any], theta: str) -> dict[str, Any]:
        branch_row = _as_mapping(self.branches.get(branch), "branch " + branch)
        cartesian = _as_mapping(branch_row.get("Cartesian"), branch + ".Cartesian")
        fields: dict[str, Any] = {field: self._reverse(cartesian.get(field), branch, field) for field in VELOCITY_FIELDS}
        pressure = cartesian.get("pressure")
        if pressure is not None:
            # The source capsule supplies pressure recipes.  Once a pressure
            # root is present it is part of this spatial jet, not a silently
            # omitted optional branch.
            fields["pressure"] = self._reverse(pressure, branch, "pressure")
        gradient: dict[str, Any] = {}
        for field in VELOCITY_FIELDS:
            physical = self._transform(fields[field]["profile_partials"], geometry, theta)
            gradient[field] = {"profile": fields[field]["profile_partials"], "physical": {axis: physical[axis] for axis in ("X", "Y", "z")}, "saved_value_ref": fields[field]["saved_value_ref"], "value_ref": fields[field]["value_ref"]}
            self.internal_counts["branch_gradient_entries"] += 3
        if "pressure" in fields:
            pressure_physical = self._transform(fields["pressure"]["profile_partials"], geometry, theta)
            gradient["pressure"] = {"profile": fields["pressure"]["profile_partials"], "physical": {axis: pressure_physical[axis] for axis in ("X", "Y", "z")}, "saved_value_ref": fields["pressure"]["saved_value_ref"], "value_ref": fields["pressure"]["value_ref"]}
            self.internal_counts["pressure_gradient_entries"] += 3
        u, v, w = gradient["u"]["physical"], gradient["v"]["physical"], gradient["w"]["physical"]
        curl = {"X": self.builder.new("add", [w["Y"], self.builder.new("negative", [v["z"]])], {"equation": "curl_X=w_Y-v_z"}), "Y": self.builder.new("add", [u["z"], self.builder.new("negative", [w["X"]])], {"equation": "curl_Y=u_z-w_X"}), "Z": self.builder.new("add", [v["X"], self.builder.new("negative", [u["Y"]])], {"equation": "curl_Z=v_X-u_Y"})}
        divergence = self.builder.new("add", [u["X"], v["Y"], w["z"]], {"equation": "divergence=u_X+v_Y+w_z", "not_numerically_zero": True, "NS_identity_not_claimed": True})
        speed = self.builder.new("add", [self.builder.new("multiply", [fields[field]["value_ref"], fields[field]["value_ref"]], {"square": field}) for field in VELOCITY_FIELDS], {"equation": "u^2+v^2+w^2"})
        vorticity = self.builder.new("add", [self.builder.new("multiply", [curl[axis], curl[axis]], {"square": "curl_" + axis}) for axis in ("X", "Y", "Z")], {"equation": "curl_X^2+curl_Y^2+curl_Z^2"})
        measure = self.builder.new("multiply", [geometry["r"]["value"], geometry["J2"]["value"]], {"equation": "r*J2", "primary_physical_measure": True})
        native_measure = self.builder.new("multiply", [measure, self._profile_R_value()], {"equation": "r*J2*R", "native_x_extra_R": True})
        integrands = {"speed_squared_physical_volume": self.builder.new("multiply", [speed, measure], {"equation": "(u^2+v^2+w^2)*(r*J2)", "not_numeric_L2": True}), "vorticity_squared_physical_volume": self.builder.new("multiply", [vorticity, measure], {"equation": "|curl|^2*(r*J2)", "not_numeric_L2": True}), "speed_squared_native_x_volume": self.builder.new("multiply", [speed, native_measure], {"equation": "(u^2+v^2+w^2)*(r*J2*R)", "native_x_extra_R": True, "not_numeric_L2": True}), "vorticity_squared_native_x_volume": self.builder.new("multiply", [vorticity, native_measure], {"equation": "|curl|^2*(r*J2*R)", "native_x_extra_R": True, "not_numeric_L2": True})}
        self.internal_counts["branch_curl_entries"] += 3
        self.internal_counts["branch_integrands"] += len(integrands)
        inverse_args = self.inverse.get("arguments", ("physical_expr:53", "physical_expr:54", "physical_expr:55", "physical_expr:56"))
        if not isinstance(inverse_args, Sequence) or isinstance(inverse_args, (str, bytes, bytearray)):
            raise KernelInputError("saved inverse arguments are malformed")
        compositions: dict[str, Any] = {}
        for field in VELOCITY_FIELDS + ("pressure",):
            item = fields.get(field)
            if _is_mapping(item) and "value_ref" in item:
                compositions[field] = self._compose(field, item["value_ref"], item["saved_value_ref"], inverse_args)
                self.internal_counts["composition_rows"] += 1
        compositions["gradient"] = {field: {axis: self._compose(field + "_" + axis, ref, fields[field]["saved_value_ref"], inverse_args) for axis, ref in gradient[field]["physical"].items()} for field in VELOCITY_FIELDS}
        if "pressure" in gradient:
            compositions["pressure_gradient"] = {axis: self._compose("pressure_" + axis, ref, fields["pressure"]["saved_value_ref"], inverse_args) for axis, ref in gradient["pressure"]["physical"].items()}
        compositions["curl"] = {axis: self._compose("curl_" + axis, ref, cartesian.get("u"), inverse_args) for axis, ref in curl.items()}
        compositions["divergence"] = self._compose("divergence", divergence, cartesian.get("u"), inverse_args)
        compositions["integrands"] = {name: self._compose(name, ref, cartesian.get("u"), inverse_args) for name, ref in integrands.items()}
        self.internal_counts["composition_rows"] += 9 + len(curl) + 1 + len(integrands) + (3 if "pressure" in gradient else 0)
        return {"branch": branch, "fields": fields, "gradient": gradient, "curl": curl, "divergence": divergence, "volume_integrands": integrands, "parametric_compositions": compositions, "geometry_shared": True, "physical_theta_not_native_phase": True, "native_fast_phase_nodes_not_differentiated": [242551, 243033], "source_values_unevaluated": True}

    def run(self) -> dict[str, Any]:
        geometry = self._geometry()
        theta = self._theta_value()
        self._emit("family_entered", {"component": COMPONENT, "new_method": "fresh reverse scalar adjoint over saved Cartesian roots"})
        frame_certificate = self._frame_zero_differential()
        self._emit("source_frame_zero_differential", {"conditional": True, "root_zero_refs": frame_certificate["root_zero_refs"], "reachable_nodes": len(frame_certificate["reachable_node_order"]), "ordinary_Z_node": frame_certificate["ordinary_Z_node"], "native_x_node": frame_certificate["native_x_node"], "fractional_part_off_seam_only": True, "N_materialized": False, "unconditional_Rm_Z_certificate": False})
        self._emit("source_join", {"saved_expression_rows": 205, "recipe_mode": self.recipe_mode, "recipe_count": len(self.recipes or {}), "fixed_port_defaults_used": not self.recipe_mode, "Rm_Z_defaulted": False, "frame_certificate_consumed": True, "source_values_unevaluated": True})
        self._emit("reverse_adjoint", {"branches": ["background", "modified"], "targets": list(VELOCITY_FIELDS), "inverse_composition_differentiated": False})
        self._emit("geometry_data", {"saved_geometry_data": True, "r": geometry["r"], "z": geometry["z"], "lambda": geometry["lambda"], "J2": geometry["J2"], "same_geometry_for_branches": True, "fresh_geometry_dual_calls": 0, "inverse": geometry["inverse"]})
        for branch in ("background", "modified"):
            if branch not in self.branches:
                raise KernelInputError("required branch missing " + branch)
            self.branch_results[branch] = self._branch(branch, geometry, theta)
        self._emit("physical_gradients", {"branches": ["background", "modified"], "gradient_entries": self.internal_counts["branch_gradient_entries"], "pressure_gradient_entries": self.internal_counts["pressure_gradient_entries"], "physical_theta_not_native_phase": True, "parametric_compositions_built": self.internal_counts["composition_rows"]})
        self._emit("curl_divergence", {"branches": ["background", "modified"], "curl_entries": self.internal_counts["branch_curl_entries"], "divergence_is_expression": True, "divergence_numeric_zero_claimed": False})
        self._emit("volume_integrands", {"branches": ["background", "modified"], "primary_measure": "r*J2", "native_x_extra_R": True, "numeric_norms_claimed": False})
        self._emit("family_completed", {"component": COMPONENT, "conditional": True, "proof_count": 1, "source_only": True, "numeric_evaluation": False})
        observed = self._event_counts()
        if observed != CALLBACK_EVENT_COUNTS:
            raise KernelInputError("callback event counts changed")
        self.internal_counts["fresh_ir_rows"] = len(self.builder.rows)
        self.internal_counts["primitive_source_ports"] = len(self.builder.primitive_cache)
        payload = {"source_saved_result": {"schema": self.saved_result.get("schema") if _is_mapping(self.saved_result) else None, "component": self.saved_result.get("component") if _is_mapping(self.saved_result) else None, "saved_pointer": "result.result.payload.expression_IR", "expression_IR_rows": 205, "not_replayed": True}, "expression_IR": list(self.builder.rows), "reverse_equation_receipts": list(self.reverse_receipts), "frame_zero_differential": frame_certificate, "primitive_source_ports": [row for row in self.builder.rows if row.get("operation") == "primitive_source_port"], "geometry": geometry, "branches": self.branch_results, "gradient_refs": {branch: {field: dict(item["physical"]) for field, item in result["gradient"].items()} for branch, result in self.branch_results.items()}, "curl_refs": {branch: result["curl"] for branch, result in self.branch_results.items()}, "divergence_refs": {branch: result["divergence"] for branch, result in self.branch_results.items()}, "integrand_refs": {branch: result["volume_integrands"] for branch, result in self.branch_results.items()}, "compositions": {branch: result["parametric_compositions"] for branch, result in self.branch_results.items()}, "composition_binding_contract": self._composition_binding_contract(), "scope": {"fixed_time": True, "profile_coordinates": ["R", "Z"], "physical_coordinates": ["X", "Y", "z"], "physical_theta": "theta", "native_fast_phase_nodes": [242551, 243033], "physical_theta_is_not_native_fast_phase": True, "primary_physical_volume_measure": "r*J2", "native_x_volume_measure": "r*J2*R", "numeric_values": False, "N_materialized": False, "Rm_Z_defaulted": False, "source_oracles": False, "global_independence": False, "axis_seam_Z4_NS": False, "source_frame_admitted_unconditionally": False, "source_frame": self.capsule.get("source_frame")}, "events": list(self.events), "source_partial_recipes_resolved_inside_kernel": self.recipe_mode, "saved_inverse_used_as_data_only": True, "source_function_bindings_recursive": True, "native_references_not_numeric_constants": True, "internal_counts": dict(self.internal_counts), "success_counts_contract": dict(SUCCESS_COUNTS), "success_flags_contract": dict(SUCCESS_FLAGS), "callback_event_counts_contract": dict(CALLBACK_EVENT_COUNTS), "count_event_map_contract": {name: list(targets) for name, targets in COUNT_EVENT_MAP.items()}}
        return {"schema": SCHEMA, "component": COMPONENT, "passed": True, "metadata_only": False, "source_only": True, "full_goal_complete": False, "numeric_evaluation": False, "runtime_model_metadata": None, "actual_counts": self._manifest_counts(observed), "new_flags": dict(SUCCESS_FLAGS), "flags": dict(SUCCESS_FLAGS), "proof_count": 1, "payload": payload, "conditional": True, "conditional_on_declared_source_atoms": True, "callback_event_count": sum(observed.values()), "callback_event_counts": observed, "failure": None}


def _rejected_report(capsule: Any, error: Exception, events: Sequence[Mapping[str, Any]] = (), internal_counts: Mapping[str, Any] | None = None) -> dict[str, Any]:
    flags = {name: False for name in SUCCESS_FLAGS}
    counts = {name: 0 for name in SUCCESS_COUNTS}
    counts["failed_family_calls"] = 1
    counts["api_refusal"] = 1
    failure = {"kind": "kernel_input_rejection", "error_type": type(error).__name__, "error": str(error), "stopped": True, "callback_prefix": list(events), "numeric_evaluation": False, "source_values_unevaluated": True}
    source_schema = capsule.get("schema") if _is_mapping(capsule) else None
    return {"schema": SCHEMA, "component": COMPONENT, "passed": False, "metadata_only": False, "source_only": True, "full_goal_complete": False, "numeric_evaluation": False, "runtime_model_metadata": None, "actual_counts": counts, "new_flags": flags, "flags": flags, "proof_count": 0, "payload": {"events": list(events), "failure": failure, "source_capsule_schema": source_schema, "internal_counts": dict(internal_counts or {}), "failure_mapping": "all_flags_false_proof_zero_manager_prefix_counts"}, "conditional": False, "failure": failure, "callback_event_count": len(events), "callback_event_counts": {event: sum(1 for item in events if item.get("event") == event) for event in CALLBACK_EVENT_COUNTS}}


def derive(capsule: Mapping[str, Any], callback: Callable[[str, Mapping[str, Any]], Any] | None = None) -> dict[str, Any]:
    """Build fresh symbolic spatial gradients, curl, and volume integrands."""

    compiler: _Compiler | None = None
    try:
        compiler = _Compiler(capsule, callback)
        return compiler.run()
    except (KernelInputError, TypeError, ValueError) as error:
        return _rejected_report(capsule, error, compiler.events if compiler is not None else (), compiler.internal_counts if compiler is not None else None)


__all__ = ["SCHEMA", "COMPONENT", "SUCCESS_COUNTS", "SUCCESS_FLAGS", "CALLBACK_EVENT_COUNTS", "COUNT_EVENT_MAP", "KernelInputError", "derive"]
