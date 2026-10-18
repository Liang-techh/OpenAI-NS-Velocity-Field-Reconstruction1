"""Bounded source-only native provenance for the actual Rm patch.

The preceding Rm native provenance result is an input receipt here.  This
module checks that receipt as data and traverses only the four new patch
source rows.  The descriptor and exact arithmetic rules come from the old
pure kernel, while the Rm controls, compact beta, and normalized gamma seed
contracts are authenticated locally before the new traversal starts.

No producer, numerical owner, accepted verifier, archive decoder, or
recurrence is called by this module.  Every returned derivative is a handle
that was already present in the authenticated source prefix.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any

try:  # The runner loads this file by a bare spec name.
    from . import original_Rm_native_derivative_provenance as _foundation
except ImportError:  # pragma: no cover - exercised by the source-only runner
    import openai_ns_reconstruction.original_Rm_native_derivative_provenance as _foundation


ProvenanceError = _foundation.ProvenanceError

DEFAULT_ORDERS = (0, 1)
_PREFIX_COUNT = 205147
_GRAPH_COUNT = 205300
_SYSTEM_HANDLE = 193470
_COORDINATE_HANDLE = 185676
_PHASE_HANDLE = 196046
_JACOBIAN_HANDLE = 187215
_INVERSE_JACOBIAN_HANDLE = 198869
_Z_HANDLE = 20
_BRANCH_RECEIPT = (
    "lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse_check.json"
)
_BRANCH_RECEIPT_SHA256 = "29781f93710407553ca958a192b326f7c9bae6ebf4f269794116c452837b0702"
_SOURCE_FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)
_COMPACT_BETA_DEFINITION = "exp(-1/(1-x^2)) for |x|<1; exact0 elsewhere"
_COMPACT_BETA_SOURCE = {
    "AST_sha256": "c4d97e103701522d471a041a3e2991031f51e5e2031d3382072cde118d1811af",
    "function": "SharedFiveMomentRepair.beta",
    "module": "lei_ren_part1_paper_compliant_five_moment_repair.py",
    "module_sha256": "978de37342cc2a63aa7c1bfda573e512f344e0157aed28b6017c7dae2ff8eba5",
}
_SOURCE_PROJECTION_BINDING = {
    "module": "lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse.py",
    "module_sha256": "b0044165763821001978135ff6601a93930a2660165d507721d23a9d9bc15a96",
    "function": "WholeZRmCenteredInverse.owner",
    "AST_sha256": "e32cb6a603c756d6cc92de64a6616c213303e3a621f0ed350e854662dba202e3",
}
_PATCH_PORTS = (
    {"field": "A", "ordinary_order": 0, "source": 197890, "target": 199051},
    {"field": "A", "ordinary_order": 1, "source": 197891, "target": 199421},
    {"field": "B", "ordinary_order": 0, "source": 197928, "target": 203304},
    {"field": "B", "ordinary_order": 1, "source": 197929, "target": 203341},
)
_LEFT_PORTS = (
    {
        "field": "A",
        "ordinary_order": 0,
        "source": 121070,
        "derivative": 178387,
        "target": 178388,
    },
    {
        "field": "A",
        "ordinary_order": 1,
        "source": 121071,
        "derivative": 178709,
        "target": 178710,
    },
    {
        "field": "B",
        "ordinary_order": 0,
        "source": 121108,
        "derivative": 182269,
        "target": 182270,
    },
    {
        "field": "B",
        "ordinary_order": 1,
        "source": 121109,
        "derivative": 182303,
        "target": 182304,
    },
)
_GAMMA_PORTS = (
    (185701, 185702, 185703),
    (185728, 185729, 185730),
    (185755, 185756, 185757),
)
_CONTROL_ROWS = (
    (193471, 193742, 193850, 193973, 194104, 194243, 194390),
    (193472, 193749, 193857, 193980, 194111, 194250, 194397),
    (193473, 193756, 193864, 193987, 194118, 194257, 194404),
    (193474, 193763, 193871, 193994, 194125, 194264, 194411),
    (193475, 193770, 193878, 194001, 194132, 194271, 194418),
)
_GRAPH_REFERENCE_KEYS = (
    "arguments",
    "operands",
    "argument",
    "numerator",
    "denominator",
    "expression",
    "variable",
    "value",
    "phase_function",
    "target_phase",
    "angle_lower",
    "angle_upper",
    "flat_predicate",
    "active_body",
    "flat_value",
    "Delta",
    "system",
    "body",
    "lower",
    "upper",
    "integrand",
    "outer",
    "dx",
    "direct",
    "chain",
    "rho",
    "Z_variable",
    "left",
    "right",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProvenanceError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    return value


def _sequence(value: Any, label: str) -> list[Any]:
    _require(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
        f"{label} must be a sequence",
    )
    return list(value)


def _handle(value: Any, label: str, count: int) -> int:
    return _foundation._handle(value, label, count)


def _row(nodes: Sequence[Any], value: Any, label: str = "node") -> Mapping[str, Any]:
    return _foundation._row(nodes, value, label)


def _op(nodes: Sequence[Any], value: int) -> str:
    return _foundation._op(nodes, value)


def _child(nodes: Sequence[Any], row: Mapping[str, Any], key: str) -> int:
    return _foundation._child(nodes, row, key)


def _normalise_digest(value: Any, label: str) -> str:
    _require(isinstance(value, str), f"{label} must be a SHA-256 digest")
    value = value.lower()
    _require(len(value) == 64, f"{label} must be a SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ProvenanceError(f"{label} must be a SHA-256 digest") from exc
    return value


def _exact(nodes: Sequence[Any], value: int) -> Any:
    return _foundation._exact(nodes, value)


def _tree_handles(value: Any, label: str, count: int) -> list[int]:
    """Collect graph handles from an explicit handle-only port structure."""

    if isinstance(value, bool):
        raise ProvenanceError(f"{label} contains a boolean handle")
    if isinstance(value, int):
        return [_handle(value, label, count)]
    if isinstance(value, Mapping):
        out: list[int] = []
        for key, child in value.items():
            out.extend(_tree_handles(child, f"{label}.{key}", count))
        return out
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        out = []
        for offset, child in enumerate(value):
            out.extend(_tree_handles(child, f"{label}[{offset}]", count))
        return out
    raise ProvenanceError(f"{label} contains a non-handle value")


def _assert_prefix_closure(
    nodes: Sequence[Mapping[str, Any]],
    roots: Sequence[int],
    prefix_count: int,
) -> None:
    """Reject selected source graphs that reach the appended history suffix."""

    seen: set[int] = set()
    stack = [int(root) for root in roots]
    while stack:
        index = stack.pop()
        _require(0 <= index < prefix_count, f"source handle {index} is outside the authenticated prefix")
        if index in seen:
            continue
        seen.add(index)
        row = _row(nodes, index)
        if row.get("operation") in (
            "exact_rational", "bound_variable", "exact_positive_integer_expression",
            "mathematical_pi", "mathematical_positive_infinity",
            "original_source_parameter", "current_original_source_parameter",
        ):
            # A rational numerator/denominator and source-parameter payloads
            # are literal data, unlike quotient and substitution child ports.
            continue
        for key in _GRAPH_REFERENCE_KEYS:
            value = row.get(key)
            if isinstance(value, int) and not isinstance(value, bool):
                stack.append(value)
            elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
                stack.extend(
                    int(child)
                    for child in value
                    if isinstance(child, int) and not isinstance(child, bool)
                )


class _PatchOracle(_foundation._Oracle):
    """Old exact lexical rules plus the authenticated patch source atoms."""

    def __init__(
        self,
        nodes: Sequence[Mapping[str, Any]],
        index: Any,
        *,
        coordinate: int,
        phase: int,
        seeds: Mapping[int, int],
        source_family: Any,
        control_nodes: Sequence[int],
        prefix_count: int,
    ) -> None:
        super().__init__(
            nodes,
            index,
            coordinate=coordinate,
            phase=phase,
            seeds=seeds,
            expected_source_family=source_family,
        )
        self.prefix_count = prefix_count
        self.control_nodes = frozenset(int(value) for value in control_nodes)
        self.control_contract_verified = False
        self.control_system_support: frozenset[str] | None = None

    def _check_compact_beta(self, index: int, row: Mapping[str, Any]) -> None:
        _require(_op(self.nodes, index) == "compact_raw_beta", f"node {index} is not compact_raw_beta")
        _require(row.get("exact_definition") == _COMPACT_BETA_DEFINITION,
                 f"compact beta {index} lacks the exact original definition")
        _require(row.get("original_radius_not_absorbed_into_raw_shape") is True,
                 f"compact beta {index} changed the original radius contract")
        _require(row.get("original_source") == _COMPACT_BETA_SOURCE,
                 f"compact beta {index} lacks the exact original-source binding")

    def _free(self, index: int) -> frozenset[str]:
        if index in self.free_cache:
            return self.free_cache[index]
        row = _row(self.nodes, index)
        operation = _op(self.nodes, index)
        if operation == "compact_raw_beta":
            self._check_compact_beta(index, row)
            result = self._free(_child(self.nodes, row, "argument"))
        elif operation == "original_Rm_unique_leading_implicit_control":
            _require(index in self.control_nodes,
                     f"unlisted unique leading control {index} is unsupported")
            _require(self.control_contract_verified,
                     "unique leading control reached before its system contract was verified")
            _require(self.control_system_support is not None,
                     "unique system input supports are not bound")
            result = self.control_system_support | self._free(_child(self.nodes, row, "Z_variable"))
        else:
            result = super()._free(index)
        self.free_cache[index] = result
        return result

    def _lookup_copied_core(self, expected: Mapping[str, Any], changed: Sequence[str]) -> int:
        candidates: list[int] = []
        for index in range(self.prefix_count):
            candidate = self.nodes[index]
            if not isinstance(candidate, Mapping) or _op(self.nodes, index) != expected.get("operation"):
                continue
            if all(
                candidate.get(key) == value
                for key, value in expected.items()
                if key not in {"operation", "op", *changed}
            ) and all(candidate.get(key) == expected.get(key) for key in changed):
                candidates.append(index)
        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise _foundation._OpenRule("copied analytic core component is absent from the authenticated prefix")
        raise _foundation._OpenRule("copied analytic core component is ambiguous in the authenticated prefix")


def _check_compact_beta_node(nodes: Sequence[Mapping[str, Any]], index: int) -> dict[str, Any]:
    row = _row(nodes, index, "compact beta")
    _require(_op(nodes, index) == "compact_raw_beta", f"gamma source {index} is not compact_raw_beta")
    _require(row.get("exact_definition") == _COMPACT_BETA_DEFINITION,
             f"compact beta {index} lacks the exact original definition")
    _require(row.get("original_radius_not_absorbed_into_raw_shape") is True,
             f"compact beta {index} changed the original radius contract")
    _require(row.get("original_source") == _COMPACT_BETA_SOURCE,
             f"compact beta {index} lacks the exact original-source binding")
    return dict(row)


def _check_gamma_lazy(
    nodes: Sequence[Mapping[str, Any]],
    index: int,
    beta: Mapping[str, Any],
) -> dict[str, Any]:
    row = _row(nodes, index, "gamma derivative")
    _require(_op(nodes, index) == "original_lazy_flat_branch",
             f"gamma derivative {index} is not an original lazy flat branch")
    _require(row.get("Delta") == beta.get("argument"),
             f"gamma derivative {index} changed the beta argument")
    _require(row.get("eta") == 1 and row.get("flat_value") == 0,
             f"gamma derivative {index} changed the original flat branch")
    _require(row.get("original_beta_exact_flat_derivative") is True,
             f"gamma derivative {index} lacks the original beta derivative contract")
    _require(row.get("active_body_not_evaluated_on_flat_branch") is True,
             f"gamma derivative {index} evaluates an original flat branch body")
    predicate = _row(nodes, _child(nodes, row, "flat_predicate"), "gamma flat predicate")
    _require(predicate.get("operation") == "logical_or" and len(predicate.get("arguments", ())) == 2,
             "gamma flat predicate is not the original two-sided cutoff")
    for ref, operator, endpoint in zip(predicate["arguments"], ("<=", ">="), (-1, 1)):
        comparison = _row(nodes, ref, "gamma flat comparison")
        _require(comparison.get("operation") == "exact_real_comparison"
                 and comparison.get("operator") == operator
                 and comparison.get("left") == beta["argument"]
                 and _exact(nodes, _child(nodes, comparison, "right")) == endpoint,
                 "gamma flat predicate changed the original support endpoints")
    return dict(row)


def _normalised_gamma_seeds(
    value: Any,
    nodes: Sequence[Mapping[str, Any]],
    prefix_count: int,
) -> tuple[dict[int, int], list[dict[str, Any]]]:
    """Authenticate the producer's normalized gamma triples without deriving them."""

    if isinstance(value, Mapping):
        value = value.get("rows", value.get("gamma", value.get("handles")))
    triples = _sequence(value, "gamma_seed_rows")
    _require(len(triples) == len(_GAMMA_PORTS), "gamma_seed_rows must contain three normalized triples")
    normalized: list[tuple[int, int, int]] = []
    for offset, triple in enumerate(triples):
        values = _sequence(triple, f"gamma_seed_rows[{offset}]")
        _require(len(values) == 3, "each gamma seed entry must contain three normalized roots")
        handles = tuple(_handle(item, f"gamma_seed_rows[{offset}][{i}]", len(nodes)) for i, item in enumerate(values))
        _require(handles == _GAMMA_PORTS[offset],
                 f"gamma seed {offset} is not the pinned normalized producer triple")
        _require(all(handle < prefix_count for handle in handles),
                 f"gamma seed {offset} is outside the authenticated prefix")
        normalized.append(handles)

    seeds: dict[int, int] = {}
    certificate: list[dict[str, Any]] = []
    scale_nodes = (2, 184121, 184123)
    scale_values = (40, 1600, 64000)
    for triple_offset, (root, first, second) in enumerate(normalized):
        for root_handle, scale_handle, scale_value in zip(
            (root, first, second), scale_nodes, scale_values
        ):
            _require(_exact(nodes, scale_handle) == scale_value,
                     f"gamma seed {triple_offset} changed its exact normalization scale")
            row = _row(nodes, root_handle, f"gamma seed {triple_offset} normalized root")
            _require(_op(nodes, root_handle) == "product",
                     f"gamma seed {triple_offset} normalized root is not a product")
            args = row.get("arguments", row.get("operands"))
            _require(args == [scale_handle, 184119, row.get("arguments", [None, None, None])[-1]],
                     f"gamma seed {triple_offset} normalized product topology changed")
        # The product arguments are intentionally read directly from the saved
        # rows, so the scale and normalization factors cannot be synthesized.
        root_args = list(_row(nodes, root, "gamma normalized beta").get("arguments", ()))
        first_args = list(_row(nodes, first, "gamma normalized first derivative").get("arguments", ()))
        second_args = list(_row(nodes, second, "gamma normalized second derivative").get("arguments", ()))
        _require(root_args[:2] == [2, 184119] and len(root_args) == 3,
                 f"gamma seed {triple_offset} beta normalization topology changed")
        _require(first_args[:2] == [184121, 184119] and len(first_args) == 3,
                 f"gamma seed {triple_offset} first normalization topology changed")
        _require(second_args[:2] == [184123, 184119] and len(second_args) == 3,
                 f"gamma seed {triple_offset} second normalization topology changed")
        norm = _row(nodes, 184119, "gamma normalization")
        _require(
            _op(nodes, 184119) == "positive_quotient"
            and norm.get("numerator") == 1
            and norm.get("denominator") == 184092
            and norm.get("source_positive_certificate") == "positive exact original beta normalization",
            "gamma normalization source contract changed",
        )
        beta = _check_compact_beta_node(nodes, root_args[2])
        first_lazy = _check_gamma_lazy(nodes, first_args[2], beta)
        second_lazy = _check_gamma_lazy(nodes, second_args[2], beta)
        # The source producer's raw beta and lazy rows are retained as data;
        # normalized roots are the actual seeds used by this patch operator.
        old_first = seeds.get(root)
        if old_first is not None and old_first != first:
            raise ProvenanceError(f"conflicting normalized gamma seed for source node {root}")
        old_second = seeds.get(first)
        if old_second is not None and old_second != second:
            raise ProvenanceError(f"conflicting normalized gamma seed for source node {first}")
        seeds[root] = first
        seeds[first] = second
        certificate.append(
            {
                "triple": [root, first, second],
                "normalization": {"beta": 2, "first": 184121, "second": 184123, "common": 184119},
                "raw_beta": root_args[2],
                "first_lazy": first_args[2],
                "second_lazy": second_args[2],
                "raw_beta_contract": {
                    "operation": beta.get("operation"),
                    "exact_definition": beta.get("exact_definition"),
                    "original_source": deepcopy(beta.get("original_source")),
                    "original_radius_not_absorbed_into_raw_shape": beta.get(
                        "original_radius_not_absorbed_into_raw_shape"
                    ),
                },
                "lazy_contract": [
                    {
                        "operation": first_lazy.get("operation"),
                        "original_beta_exact_flat_derivative": first_lazy.get(
                            "original_beta_exact_flat_derivative"
                        ),
                        "active_body_not_evaluated_on_flat_branch": first_lazy.get(
                            "active_body_not_evaluated_on_flat_branch"
                        ),
                    },
                    {
                        "operation": second_lazy.get("operation"),
                        "original_beta_exact_flat_derivative": second_lazy.get(
                            "original_beta_exact_flat_derivative"
                        ),
                        "active_body_not_evaluated_on_flat_branch": second_lazy.get(
                            "active_body_not_evaluated_on_flat_branch"
                        ),
                    },
                ],
            }
        )
    return seeds, certificate


def _validate_system(
    nodes: Sequence[Mapping[str, Any]],
    index: Any,
    leading: Mapping[str, Any],
    source_family: Any,
    branch_receipt_sha256: str,
    oracle: _PatchOracle,
) -> tuple[list[int], dict[str, Any]]:
    """Bind the original five-function system and prove its lexical supports."""

    _require(leading.get("system") == _SYSTEM_HANDLE,
             "leading_bindings.system must be the authenticated original system 193470")
    system = _row(nodes, leading["system"], "original leading implicit system")
    _require(_op(nodes, _SYSTEM_HANDLE) == "original_Rm_unique_leading_implicit_system",
             "leading system operation changed")
    _require(system.get("source_family") == source_family, "leading system source family changed")
    if leading.get("source_family") is not None:
        _require(leading.get("source_family") == source_family, "leading binding source family changed")
    _require(system.get("defining_map") == "L*h+Q(h,invAm2)+true_Rm_defects=0",
             "leading system defining map changed")
    _require(system.get("branch_receipt_sha256") == branch_receipt_sha256 == _BRANCH_RECEIPT_SHA256,
             "leading system branch receipt is not the pinned original receipt")
    _require(system.get("original_common_box_branch_receipt") == _BRANCH_RECEIPT,
             "leading system branch receipt filename changed")
    _require(system.get("original_whole_Z_base_unique_branch_and_Jacobian_admitted") is True,
             "original unique branch/Jacobian admission is missing")
    _require(system.get("original_leading_fixed_point_not_Rc_finite_N_control") is True,
             "leading fixed point was relabeled as an Rm finite-N control")
    _require(system.get("numerical_enclosures_not_values") is True,
             "leading system was changed into numeric values")
    _require(system.get("source_projection_binding") == _SOURCE_PROJECTION_BINDING,
             "leading system source projection binding changed")
    _require(system.get("source_bindings") == leading.get("bindings"),
             "leading system source bindings changed")
    branch_transfer = leading.get("branch_transfer")
    _require(isinstance(branch_transfer, Mapping), "leading branch transfer metadata is required")
    _require(system.get("exact_source_and_weight_branch_transfer") == branch_transfer,
             "leading system no longer carries the admitted branch transfer")
    _require(branch_transfer.get("original_whole_Z_unique_branch_receipt") == _BRANCH_RECEIPT,
             "leading branch transfer receipt changed")
    _require(branch_transfer.get("actual_leading_N0_equal_by_exact_integer") is True,
             "leading branch transfer lost the fixed-input identity")
    _require(branch_transfer.get("exact_defects_from_original_physical_primitives") is True,
             "leading branch transfer lost exact defect provenance")
    _require(branch_transfer.get("original_raw_beta_definition_and_finite_support_normalization_unchanged") is True,
             "leading branch transfer lost raw beta provenance")
    _require(branch_transfer.get("original_four_cell_C0_and_Jacobian_contractions_strictly_lt1") is True,
             "leading branch transfer lost its Jacobian admission")

    weights = leading.get("W")
    defects = leading.get("defects")
    alpha = leading.get("alpha")
    units = leading.get("units")
    _require(system.get("exact_weight_nodes") == weights, "leading system exact weights changed")
    _require(
        system.get("exact_defect_value_nodes") == [row[0] for row in _sequence(defects, "leading defects")],
        "leading system exact defect roots changed",
    )
    _require(system.get("invAm2") == _sequence(alpha, "leading alpha")[0],
             "leading system exact invAm2 changed")
    _require(system.get("exact_Z_independent_row_and_control_units") == units,
             "leading system exact units changed")

    controls = _sequence(leading.get("controls"), "leading controls")
    _require(len(controls) == 5, "leading controls must have five components")
    control_nodes: list[int] = []
    for component, raw_row in enumerate(controls):
        row_handles = _sequence(raw_row, f"leading controls[{component}]")
        _require(len(row_handles) == 7, "each leading control row must contain seven ordinary entries")
        handles = [_handle(item, f"leading controls[{component}][{order}]", len(nodes)) for order, item in enumerate(row_handles)]
        _require(tuple(handles) == _CONTROL_ROWS[component], "actual higher control source ports changed")
        _require(all(item < _PREFIX_COUNT for item in handles),
                 f"leading controls[{component}] reaches outside the authenticated prefix")
        value_handle = handles[0]
        _require(value_handle == 193471 + component,
                 f"leading control component {component} value root changed")
        value_row = _row(nodes, value_handle, f"leading control value {component}")
        _require(_op(nodes, value_handle) == "original_Rm_unique_leading_implicit_control",
                 f"leading control value {component} operation changed")
        _require(
            value_row.get("system") == _SYSTEM_HANDLE
            and value_row.get("component") == component
            and value_row.get("Z_variable") == _Z_HANDLE
            and value_row.get("source_family") == source_family
            and value_row.get("exact_value_of_same_unique_original_common_box_solution") is True,
            f"leading control value {component} contract changed",
        )
        control_nodes.append(value_handle)

    _require(_op(nodes, _Z_HANDLE) == "bound_variable" and _row(nodes, _Z_HANDLE).get("name") == "Z",
             "leading controls must use the original Z variable")
    _require(_tree_handles(weights, "leading W", _PREFIX_COUNT), "leading W is empty")
    _require(len(_sequence(weights.get("L"), "leading W.L")) == 5, "leading W.L must have five rows")
    for key, length in (("fg", 2), ("gg", 2), ("ff", 3), ("ff_over_x", 3)):
        _require(len(_sequence(weights.get(key), f"leading W.{key}")) == length,
                 f"leading W.{key} shape changed")
    L = _sequence(weights.get("L"), "leading W.L")
    _require(all(len(_sequence(row, "leading W.L row")) == 5 for row in L), "leading W.L shape changed")
    _require(len(_sequence(defects, "leading defects")) == 5, "leading defects must have five rows")
    _require(all(len(_sequence(row, "leading defect row")) == 7 for row in defects),
             "leading defect rows must contain seven ordinary entries")
    _require(len(_sequence(alpha, "leading alpha")) == 7, "leading alpha must contain seven entries")
    _require(len(_sequence(units, "leading units")) == 5, "leading units must contain five entries")

    weight_supports: list[dict[str, Any]] = []
    for key, value in weights.items():
        for handle in _tree_handles(value, f"leading W.{key}", _PREFIX_COUNT):
            free = oracle._free(handle)
            _require(not free, f"leading weight {handle} depends on a free variable: {sorted(free)}")
            weight_supports.append({"path": key, "node": handle, "free_variables": []})
    defect_supports: list[dict[str, Any]] = []
    for component, row in enumerate(defects):
        for order, item in enumerate(row):
            handle = _handle(item, f"leading defects[{component}][{order}]", _PREFIX_COUNT)
            free = oracle._free(handle)
            _require(free <= frozenset({"Z"}),
                     f"leading defect {handle} has unsupported free variables: {sorted(free)}")
            defect_supports.append({"component": component, "ordinary_order": order, "node": handle,
                                    "free_variables": sorted(free)})
    alpha_supports: list[dict[str, Any]] = []
    for order, item in enumerate(alpha):
        handle = _handle(item, f"leading alpha[{order}]", _PREFIX_COUNT)
        free = oracle._free(handle)
        _require(free <= frozenset({"Z"}),
                 f"leading alpha {handle} has unsupported free variables: {sorted(free)}")
        alpha_supports.append({"ordinary_order": order, "node": handle, "free_variables": sorted(free)})
    unit_supports: list[dict[str, Any]] = []
    for order, item in enumerate(units):
        handle = _handle(item, f"leading units[{order}]", _PREFIX_COUNT)
        free = oracle._free(handle)
        _require(not free, f"leading unit {handle} depends on a free variable: {sorted(free)}")
        unit_supports.append({"ordinary_order": order, "node": handle, "free_variables": []})

    oracle.control_nodes = frozenset(control_nodes)
    system_support = frozenset().union(*(frozenset(row["free_variables"])
        for row in weight_supports + defect_supports + alpha_supports + unit_supports))
    oracle.control_system_support = system_support
    oracle.control_contract_verified = True
    control_supports: list[dict[str, Any]] = []
    for component, raw_row in enumerate(controls):
        for order, item in enumerate(raw_row):
            handle = _handle(item, f"leading controls[{component}][{order}]", _PREFIX_COUNT)
            free = oracle._free(handle)
            _require(free <= frozenset({"Z"}),
                     f"leading control row {component},{order} has unsupported free variables: {sorted(free)}")
            _require(order != 0 or free == frozenset({"Z"}),
                     f"leading control value {component} lost explicit Z dependence")
            control_supports.append({
                "component": component,
                "ordinary_order": order,
                "node": handle,
                "free_variables": sorted(free),
                "native_coordinate_independent_by_lexical_support": "true_Rm_patch_x" not in free,
            })

    certificate = {
        "passed": True,
        "system": _SYSTEM_HANDLE,
        "source_family": deepcopy(source_family),
        "branch_receipt_sha256": branch_receipt_sha256,
        "control_value_nodes": list(control_nodes),
        "control_support_rows": control_supports,
        "weight_supports": weight_supports,
        "defect_supports": defect_supports,
        "alpha_supports": alpha_supports,
        "unit_supports": unit_supports,
        "all_higher_control_rows_checked_by_lexical_dependency": True,
        "system_input_free_variables": sorted(system_support),
        "implicit_recurrence_replayed": False,
        "unique_branch_support_theorem": "same explicit system inputs and Z give the same unique control; no native parameter is present",
        "control_values_not_set_to_zero": True,
        "same_unique_original_common_box_solution": True,
        "Z_dependence_retained": True,
    }
    return control_nodes, certificate


def _validate_accepted_left(
    value: Any,
    expected_nodes_digest: str,
    derivative_prefix_digest: str,
) -> Mapping[str, Any]:
    left = _mapping(value, "accepted_left_provenance")
    _require(left.get("passed") is True and left.get("status") == "proved",
             "accepted left provenance must contain four passed rows")
    for key in (
        "source_only",
        "source_executed",
        "actual_saved_graph_run",
        "synthetic_execution",
        "nodes_input_unmutated",
        "closure_claim",
        "numeric_values_evaluated",
    ):
        expected = key in {"source_only", "nodes_input_unmutated"}
        _require(left.get(key) is expected, f"accepted left provenance flag {key} changed")
    _require(left.get("side") == "postslow" and left.get("ordinary_orders") == [0, 1],
             "accepted left provenance scope changed")
    _require(left.get("expected_nodes_digest") == expected_nodes_digest == left.get("nodes_digest"),
             "accepted left graph digest is not bound to the new graph")
    _require(left.get("derivative_source_prefix_count") == _PREFIX_COUNT,
             "accepted left derivative prefix count changed")
    _require(left.get("derivative_source_prefix_digest") == derivative_prefix_digest,
             "accepted left derivative prefix digest changed")
    _require(left.get("coordinate_handle") == 84528 and left.get("actual_phase_handle") == 119416,
             "accepted left source ports changed")
    _require(left.get("physical_J_handle") == 1 and left.get("inverse_J_handle") == 172066,
             "accepted left Jacobian ports changed")
    _require(left.get("unresolved_gates") == [] and left.get("closure_claim") is False,
             "accepted left unresolved gates changed")
    rows = _sequence(left.get("rows"), "accepted left rows")
    _require(len(rows) == len(_LEFT_PORTS), "accepted left provenance must contain exactly four rows")
    for expected, row in zip(_LEFT_PORTS, rows):
        row = _mapping(row, "accepted left row")
        _require(
            row.get("side") == "postslow"
            and row.get("field") == expected["field"]
            and row.get("ordinary_order") == expected["ordinary_order"]
            and row.get("source_handle") == expected["source"]
            and row.get("source_node") == expected["source"]
            and row.get("physical_target_handle") == expected["target"]
            and row.get("target_handle") == expected["target"]
            and row.get("derivative_handle") == expected["derivative"]
            and row.get("derivative_root") == expected["derivative"]
            and row.get("expected_physical_target_handle") == expected["target"]
            and row.get("physical_target_matches") is True
            and row.get("status") == "proved",
            f"accepted left source/target port changed for {expected['field']} Z{expected['ordinary_order']}",
        )
    visited = left.get("visited_source_D_pairs")
    _require(isinstance(visited, Sequence) and not isinstance(visited, (str, bytes, bytearray)),
             "accepted left visited source pairs must be retained as data")
    _require(len(visited) == 365, "accepted left visited source-pair receipt changed")
    _require(isinstance(left.get("visited_source_to_D_pairs"), Mapping),
             "accepted left source-to-derivative map must be retained as data")
    return left


def _source_seed_handles(
    nodes: Sequence[Mapping[str, Any]],
    leading: Mapping[str, Any],
    radial_background: Mapping[str, Any],
    coordinate: int,
    phase: int,
    jacobian: int,
    index: Any,
    source_family: Any,
    gamma_seed_rows: Any,
    prefix_count: int,
) -> tuple[dict[int, int], dict[str, Any], list[dict[str, Any]]]:
    chart = _mapping(leading.get("chart"), "leading_bindings.chart")
    _require(chart.get("phase") == coordinate,
             "leading chart native coordinate does not match the patch coordinate")
    _require(chart.get("physical_measure") == jacobian,
             "leading chart Jacobian does not match the patch chart")
    leading_handles, leading_descriptors = _foundation._source_rows(chart, radial_background, nodes)
    _require(leading_handles and leading_descriptors, "leading source field germs are missing")
    _require(all(order in (0, 1, 2) for _, order in leading_handles),
             "leading source germs must be restricted to ordinary orders 0..2")

    gamma_seeds, gamma_certificate = _normalised_gamma_seeds(gamma_seed_rows, nodes, prefix_count)
    seed_oracle = _PatchOracle(
        nodes,
        index,
        coordinate=coordinate,
        phase=phase,
        seeds={},
        source_family=source_family,
        control_nodes=(),
        prefix_count=prefix_count,
    )
    zero = index.rational(0)
    seeds: dict[int, int] = {phase: zero}
    for (field, order), source_handle in leading_handles.items():
        radial_handle = leading_descriptors[field]["radial"][order]
        derivative_handle = seed_oracle._mul(jacobian, radial_handle)
        _require(source_handle not in (zero, index.rational(1)) or derivative_handle == zero,
                 f"constant leading field seed {source_handle} has a nonzero native derivative")
        previous = seeds.get(source_handle)
        if previous is not None and previous != derivative_handle:
            raise ProvenanceError(f"conflicting leading-field seed for source node {source_handle}")
        seeds[source_handle] = derivative_handle
    for source_handle, derivative_handle in gamma_seeds.items():
        _require(source_handle not in (zero, index.rational(1)) or derivative_handle == zero,
                 f"constant normalized gamma seed {source_handle} has a nonzero native derivative")
        previous = seeds.get(source_handle)
        if previous is not None and previous != derivative_handle:
            raise ProvenanceError(f"conflicting normalized gamma seed for source node {source_handle}")
        seeds[source_handle] = derivative_handle
    return seeds, {
        "leading_field_orders": [0, 1, 2],
        "source_family": deepcopy(source_family),
        "chart_coordinate": coordinate,
        "actual_phase": phase,
        "jacobian": jacobian,
        "leading_fields": {
            field: {
                "source": {str(order): handle for order, handle in sorted(values["source"].items())},
                "radial": {str(order): handle for order, handle in sorted(values["radial"].items())},
            }
            for field, values in leading_descriptors.items()
        },
        "gamma": gamma_certificate,
    }, gamma_certificate


def _visited(oracle: _PatchOracle) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    pairs = [
        {
            "source_handle": source,
            "derivative_context": list(context),
            "derivative_handle": derivative,
        }
        for (source, context), derivative in sorted(
            oracle.visited.items(), key=lambda item: (item[0][0], repr(item[0][1]))
        )
    ]
    mapping: dict[str, list[dict[str, Any]]] = {}
    for pair in pairs:
        mapping.setdefault(str(pair["source_handle"]), []).append(pair)
    return pairs, mapping


def derive(
    nodes: Sequence[Mapping[str, Any]],
    leading_bindings: Mapping[str, Any],
    loop_chart: Mapping[str, Any],
    radial_background: Mapping[str, Any],
    accepted_left_provenance: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    derivative_prefix_digest: str,
    source_family: Any,
    gamma_seed_rows: Sequence[Any],
    branch_receipt_sha256: str,
) -> dict[str, Any]:
    """Derive the four actual patch slow rows in the authenticated source DAG."""

    base = _sequence(nodes, "nodes")
    _require(len(base) == _GRAPH_COUNT, "current source graph must contain exactly 205300 nodes")
    _require(all(isinstance(row, Mapping) for row in base), "nodes must contain mappings")
    leading = _mapping(leading_bindings, "leading_bindings")
    loop = _mapping(loop_chart, "loop_chart")
    background = _mapping(radial_background, "radial_background")
    expected = _normalise_digest(expected_nodes_digest, "expected_nodes_digest")
    prefix_expected = _normalise_digest(derivative_prefix_digest, "derivative_prefix_digest")
    branch_digest = _normalise_digest(branch_receipt_sha256, "branch_receipt_sha256")
    before = _foundation.canonical_digest(base)
    _require(before == expected, "nodes digest does not match expected_nodes_digest")
    lookup_nodes = base[:_PREFIX_COUNT]
    lookup_digest = _foundation.canonical_digest(lookup_nodes)
    _require(lookup_digest == prefix_expected, "source prefix digest is not authenticated")
    index = _foundation._DescriptorIndex(lookup_nodes)

    accepted_left = _validate_accepted_left(
        accepted_left_provenance,
        expected,
        prefix_expected,
    )
    coordinate, phase, loop_sources, radial_targets = _foundation._loop_rows(loop, base, DEFAULT_ORDERS)
    _require(coordinate == _COORDINATE_HANDLE and phase == _PHASE_HANDLE,
             "patch coordinate or actual phase port changed")
    _require(loop.get("physical_measure") == _JACOBIAN_HANDLE
             and loop.get("invJ") == _INVERSE_JACOBIAN_HANDLE,
             "patch chart Jacobian ports changed")
    radial = _mapping(loop.get("radial"), "loop_chart.radial")
    _require(radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True,
             "patch radial source did not retain exact native/Jacobian normalization")
    jacobian, inverse_jacobian = _foundation._find_jacobian(
        base,
        index,
        _mapping(leading.get("chart"), "leading_bindings.chart"),
        loop,
    )
    _require(jacobian == _JACOBIAN_HANDLE and inverse_jacobian == _INVERSE_JACOBIAN_HANDLE,
             "patch Jacobian is not the authenticated exact 1/x pair")
    _require(loop_sources["A"][0] == 197890 and loop_sources["A"][1] == 197891
             and loop_sources["B"][0] == 197928 and loop_sources["B"][1] == 197929,
             "patch source ports changed")
    _require(radial_targets["A"][0] == 199051 and radial_targets["A"][1] == 199421
             and radial_targets["B"][0] == 203304 and radial_targets["B"][1] == 203341,
             "patch target ports changed")

    # The control oracle is deliberately unable to expose control support
    # until all independent system-input supports have been checked.
    provisional = _PatchOracle(
        base,
        index,
        coordinate=coordinate,
        phase=phase,
        seeds={},
        source_family=source_family,
        control_nodes=(),
        prefix_count=_PREFIX_COUNT,
    )
    # Validate the system first to obtain the five authenticated control roots.
    # The seed oracle is assembled immediately afterwards from the source-only
    # field germs and normalized gamma producer contract.
    control_nodes, control_certificate = _validate_system(
        base,
        index,
        leading,
        source_family,
        branch_digest,
        provisional,
    )
    seeds, source_seed_bindings, gamma_certificate = _source_seed_handles(
        base,
        leading,
        background,
        coordinate,
        phase,
        jacobian,
        index,
        source_family,
        gamma_seed_rows,
        _PREFIX_COUNT,
    )
    # Re-run only the source-input support check against the final oracle's
    # authenticated control set; no accepted left row is an oracle root.
    oracle = _PatchOracle(
        base,
        index,
        coordinate=coordinate,
        phase=phase,
        seeds=seeds,
        source_family=source_family,
        control_nodes=control_nodes,
        prefix_count=_PREFIX_COUNT,
    )
    oracle.control_contract_verified = True
    oracle.control_system_support = frozenset(control_certificate["system_input_free_variables"])

    closure_roots: list[int] = [coordinate, phase, jacobian, inverse_jacobian]
    closure_roots.extend(loop_sources[field][order] for field in ("A", "B") for order in DEFAULT_ORDERS)
    closure_roots.extend(radial_targets[field][order] for field in ("A", "B") for order in DEFAULT_ORDERS)
    closure_roots.extend(seeds)
    closure_roots.extend(seeds.values())
    closure_roots.extend(control_nodes)
    closure_roots.extend(_tree_handles(leading.get("controls"), "leading controls", _PREFIX_COUNT))
    closure_roots.extend(_tree_handles(leading.get("W"), "leading W", _PREFIX_COUNT))
    closure_roots.extend(_tree_handles(leading.get("defects"), "leading defects", _PREFIX_COUNT))
    closure_roots.extend(_tree_handles(leading.get("alpha"), "leading alpha", _PREFIX_COUNT))
    closure_roots.extend(_tree_handles(leading.get("units"), "leading units", _PREFIX_COUNT))
    for gamma_row in gamma_certificate:
        closure_roots.extend(
            gamma_row[key]
            for key in ("raw_beta", "first_lazy", "second_lazy")
            if isinstance(gamma_row.get(key), int)
        )
    _assert_prefix_closure(base, closure_roots, _PREFIX_COUNT)

    rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for port in _PATCH_PORTS:
        source_handle = port["source"]
        target_handle = port["target"]
        row: dict[str, Any] = {
            "side": "patch",
            "field": port["field"],
            "ordinary_order": port["ordinary_order"],
            "source_handle": source_handle,
            "source_node": source_handle,
            "physical_target_handle": target_handle,
            "target_handle": target_handle,
            "status": "unresolved",
            "formula": "physical_target = invJ * D_native(source) at held actual phase",
        }
        try:
            derivative_handle = oracle.derivative(source_handle)
            expected_target = oracle._mul(inverse_jacobian, derivative_handle)
            row.update(
                derivative_handle=derivative_handle,
                derivative_root=derivative_handle,
                expected_physical_target_handle=expected_target,
                physical_target_matches=target_handle == expected_target,
            )
            if target_handle != expected_target:
                row["error"] = "saved patch slow target is not invJ*D_native(source)"
                unresolved.append(
                    {
                        "gate": "physical_patch_target_descriptor_mismatch",
                        "field": port["field"],
                        "ordinary_order": port["ordinary_order"],
                        "source_handle": source_handle,
                        "physical_target_handle": target_handle,
                        "expected_physical_target_handle": expected_target,
                    }
                )
            else:
                row["status"] = "proved"
        except (ProvenanceError, _foundation._OpenRule, KeyError, TypeError, ValueError) as exc:
            row["error"] = str(exc)
            unresolved.append(
                {
                    "gate": "patch_native_derivative_rule_unresolved",
                    "field": port["field"],
                    "ordinary_order": port["ordinary_order"],
                    "source_handle": source_handle,
                    "error": str(exc),
                }
            )
        rows.append(row)

    after = _foundation.canonical_digest(base)
    if after != before:
        unresolved.append({"gate": "nodes_input_mutated", "before": before, "after": after})
    visited_pairs, visited_map = _visited(oracle)
    control_certificate = deepcopy(control_certificate)
    control_certificate["passed"] = not any(
        row.get("free_variables") and not set(row["free_variables"]) <= {"Z"}
        for row in control_certificate.get("control_support_rows", [])
    )
    control_certificate["gamma_seed_rows"] = deepcopy(gamma_certificate)
    passed = bool(rows) and not unresolved and control_certificate["passed"] and all(
        row.get("status") == "proved" for row in rows
    )
    return {
        "schema": "Rm_control_bound_patch_native_provenance_v1",
        "passed": passed,
        "status": "proved" if passed else "unresolved",
        "source_only": True,
        "source_executed": False,
        "actual_saved_graph_run": False,
        "synthetic_execution": False,
        "nodes_input_unmutated": after == before,
        "expected_nodes_digest": expected,
        "nodes_digest": after,
        "derivative_source_prefix_count": _PREFIX_COUNT,
        "derivative_source_prefix_digest": lookup_digest,
        "side": "patch",
        "ordinary_orders": list(DEFAULT_ORDERS),
        "coordinate_handle": coordinate,
        "actual_phase_handle": phase,
        "physical_J_handle": jacobian,
        "inverse_J_handle": inverse_jacobian,
        "seed_handles": {str(source): derivative for source, derivative in sorted(seeds.items())},
        "rows": rows,
        "visited_source_D_pairs": visited_pairs,
        "visited_source_to_D_pairs": visited_map,
        "unresolved_gates": unresolved,
        "closure_claim": False,
        "cross_chart_native_seam_proved": False,
        "cross_chart_seam_claim": False,
        "numeric_values_evaluated": False,
        "numeric_field_evaluated": False,
        "accepted_left_provenance": deepcopy(accepted_left),
        "accepted_left_provenance_digest": _foundation.canonical_digest(accepted_left),
        "accepted_left_provenance_consumed_as_data": True,
        "accepted_left_traversal_replayed": False,
        "accepted_left_visited_source_D_pairs": deepcopy(accepted_left["visited_source_D_pairs"]),
        "accepted_left_visited_source_to_D_pairs": deepcopy(accepted_left["visited_source_to_D_pairs"]),
        "accepted_left_unresolved_gates": deepcopy(accepted_left["unresolved_gates"]),
        "source_bindings": deepcopy(leading.get("bindings")),
        "source_seed_bindings": deepcopy(source_seed_bindings),
        "control_certificate": control_certificate,
        "control_support_rows": deepcopy(control_certificate["control_support_rows"]),
        "free_supports": {
            "weights": deepcopy(control_certificate["weight_supports"]),
            "defects": deepcopy(control_certificate["defect_supports"]),
            "alpha": deepcopy(control_certificate["alpha_supports"]),
            "units": deepcopy(control_certificate["unit_supports"]),
            "controls": deepcopy(control_certificate["control_support_rows"]),
        },
        "gamma_seed_certificate": deepcopy(gamma_certificate),
        "branch_receipt_sha256": branch_digest,
        "source_family": deepcopy(source_family),
        "no_added_nodes": True,
    }


__all__ = ["DEFAULT_ORDERS", "ProvenanceError", "derive"]
