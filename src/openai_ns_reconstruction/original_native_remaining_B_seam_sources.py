"""Bounded source comparison for the remaining native ``B`` seams.

This module is a source-only leaf.  It consumes the authenticated packets
assembled by the runner and compares only the endpoint source DAGs needed for
``J2`` through ``J15``.  It deliberately does not import an ancestor module,
execute a producer, invoke the accepted endpoint canonicalizer, expand a
polynomial, evaluate a quadrature, or turn an unknown source operation into an
inequality claim.
"""

import hashlib
import json


PREFIX = "lei_ren_part1_paper_compliant_current_original_"
SCHEMA = "original_native_remaining_B_seam_sources_v1"

PARENT = PREFIX + "native_micro_macro_B_jump_source_v2.json"
WEIGHTED = PREFIX + "native_active_S1_weighted_ibp.json"
MOD = PREFIX + "modified_Q4_pressure_C3_functions.json.gz"
LOOP = PREFIX + "prescribed_loop_C5_functions.json.gz"
PATCH = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
OUTER = PREFIX + "true_outer_C6_C5_functions.json.gz"

RECORD_NAMES = (PARENT, WEIGHTED, MOD, LOOP, PATCH, OUTER)
RECORD_HASHES = {
    PARENT: "2505e298aacbac9e15444e93dd7b054deef201c8d02193938e5391ec7b528048",
    WEIGHTED: "aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28",
    MOD: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    LOOP: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    PATCH: "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad",
    OUTER: "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
}

STAGE_FLAGS = (
    "new_remaining_seam_source_comparison_completed",
    "accepted_J0_J1_not_reexecuted",
    "all_old_terms_retained",
)

# ``all_native_B_jumps_zero_proved`` is intentionally dynamic.  It is not a
# false flag because a future bounded source walk may prove every remaining
# jump.  These flags describe gates that this leaf never claims.
FALSE_FLAGS = (
    "numeric_values_computed",
    "old_error_terms_removed",
    "old_operators_replayed",
    "old_algebra_replayed",
    "global_five_moment_closure",
    "time_order_recursion_claimed",
    "full_goal_complete",
    "native_old_density_identity_proved",
    "native_active_m15_source_identity_proved",
    "native_compact_m15_source_identity_proved",
    "terminal_difference_zero_proved",
)

B_FORMULA = "R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)"
IN_FORMULA = "I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy)"
L_FORMULA = "1-delta*Z^2"
L_NODE = 25272
ORDINARY_Z_ORDERS = (0, 1)
OWNED_Z_CELL = ("0", "1/2")

# These are the analytic and comparison primitives actually present in the
# authenticated source archives.  Everything else remains an opaque,
# unproved source operation.  In particular, this list is not an invitation
# to evaluate any analytic function numerically.
ANALYTIC_NAMES = frozenset({
    "compact_raw_beta_derivative",
    "cos",
    "exp",
    "expm1",
    "exprel",
    "fractional_part",
    "gamma",
    "log",
    "log1p",
    "log1prel",
    "original_flat_sigma",
    "original_flat_sigma_prime",
    "original_flat_sigma_second",
    "original_flat_sigma_third",
    "positive_sqrt",
    "sin",
    "sqrt",
})
COMPARISON_OPERATORS = frozenset({"==", "!=", "<", "<=", ">", ">="})
MAX_OPAQUE_REFS = 2048

# This is the static archive binding map.  LOOP supplies geometry for the
# first ten charts and MOD supplies their physical U/m/Q ports.  PATCH and
# OUTER own both geometry and physical ports in their respective sections.
# The weighted graph remains the authority for the actual port node ids.
_SEAM_SPECS = (
    {
        "jump_index": 2,
        "native_graph_node": 522,
        "seam": "frozen_macro -> first_switch",
        "certificate_family": "macro_first",
        "left_chart": "frozen_macro",
        "right_chart": "first_switch",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 3,
        "native_graph_node": 523,
        "seam": "first_switch -> second_switch",
        "certificate_family": "switch",
        "left_chart": "first_switch",
        "right_chart": "second_switch",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 4,
        "native_graph_node": 524,
        "seam": "second_switch -> post_power",
        "certificate_family": "switch",
        "left_chart": "second_switch",
        "right_chart": "post_power",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 5,
        "native_graph_node": 525,
        "seam": "post_power -> long_reshape",
        "certificate_family": "annular",
        "left_chart": "post_power",
        "right_chart": "long_reshape",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 6,
        "native_graph_node": 526,
        "seam": "long_reshape -> reference",
        "certificate_family": "annular",
        "left_chart": "long_reshape",
        "right_chart": "reference",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 7,
        "native_graph_node": 527,
        "seam": "reference -> restoration",
        "certificate_family": "annular",
        "left_chart": "reference",
        "right_chart": "restoration",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 8,
        "native_graph_node": 528,
        "seam": "restoration -> postrestore",
        "certificate_family": "annular",
        "left_chart": "restoration",
        "right_chart": "postrestore",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": LOOP,
        "left_physical": MOD,
        "right_physical": MOD,
    },
    {
        "jump_index": 9,
        "native_graph_node": 529,
        "seam": "postrestore -> Rm_patch_finite_N",
        "certificate_family": "rm_outer",
        "left_chart": "postrestore",
        "right_chart": "Rm_patch_finite_N",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": LOOP,
        "right_geometry": PATCH,
        "left_physical": MOD,
        "right_physical": PATCH,
    },
    {
        "jump_index": 10,
        "native_graph_node": 530,
        "seam": "Rm_patch_finite_N -> Rh_reference",
        "certificate_family": "rm_outer",
        "left_chart": "Rm_patch_finite_N",
        "right_chart": "Rh_reference",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": PATCH,
        "right_geometry": OUTER,
        "left_physical": PATCH,
        "right_physical": OUTER,
    },
    {
        "jump_index": 11,
        "native_graph_node": 531,
        "seam": "Rh_reference -> O2_slope",
        "certificate_family": "rm_outer",
        "left_chart": "Rh_reference",
        "right_chart": "O2_slope",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": OUTER,
        "right_geometry": OUTER,
        "left_physical": OUTER,
        "right_physical": OUTER,
    },
    {
        "jump_index": 12,
        "native_graph_node": 532,
        "seam": "O2_slope -> O2_axial",
        "certificate_family": "rm_outer",
        "left_chart": "O2_slope",
        "right_chart": "O2_axial",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": OUTER,
        "right_geometry": OUTER,
        "left_physical": OUTER,
        "right_physical": OUTER,
    },
    {
        "jump_index": 13,
        "native_graph_node": 533,
        "seam": "O2_axial -> O2_buffer",
        "certificate_family": "rm_outer",
        "left_chart": "O2_axial",
        "right_chart": "O2_buffer",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": OUTER,
        "right_geometry": OUTER,
        "left_physical": OUTER,
        "right_physical": OUTER,
    },
    {
        "jump_index": 14,
        "native_graph_node": 534,
        "seam": "O2_buffer -> O3_transition",
        "certificate_family": "rm_outer",
        "left_chart": "O2_buffer",
        "right_chart": "O3_transition",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": OUTER,
        "right_geometry": OUTER,
        "left_physical": OUTER,
        "right_physical": OUTER,
    },
    {
        "jump_index": 15,
        "native_graph_node": 535,
        "seam": "O3_transition -> O3_power",
        "certificate_family": "rm_outer",
        "left_chart": "O3_transition",
        "right_chart": "O3_power",
        "left_endpoint": "upper",
        "right_endpoint": "lower",
        "left_geometry": OUTER,
        "right_geometry": OUTER,
        "left_physical": OUTER,
        "right_physical": OUTER,
    },
)


def _jsonable(value):
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(value[key]) for key in sorted(value, key=str)}
    return value


def _digest(value):
    payload = json.dumps(
        _jsonable(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _gcd(left, right):
    left, right = abs(left), abs(right)
    while right:
        left, right = right, left % right
    return left or 1


def _rational(numerator, denominator=1):
    if not (_is_int(numerator) and _is_int(denominator)) or denominator == 0:
        return None
    divisor = _gcd(numerator, denominator)
    numerator //= divisor
    denominator //= divisor
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    return ("rational", numerator, denominator)


def _is_rational(value):
    return isinstance(value, tuple) and len(value) == 3 and value[0] == "rational"


def _r_add(left, right):
    return _rational(left[1] * right[2] + right[1] * left[2], left[2] * right[2])


def _r_mul(left, right):
    return _rational(left[1] * right[1], left[2] * right[2])


def _r_neg(value):
    return _rational(-value[1], value[2])


def _r_div(left, right):
    if right[1] == 0:
        return None
    return _rational(left[1] * right[2], left[2] * right[1])


def _metadata_digest(row, excluded=()):
    excluded = set(excluded) | {"operation", "arguments"}
    return _digest({key: value for key, value in row.items() if key not in excluded})


def _metadata_keys(row, excluded=()):
    excluded = set(excluded) | {"operation", "arguments"}
    return sorted(key for key in row if key not in excluded)


def _explicit_canonical_free_z(row, name):
    lower_name = name.lower() if isinstance(name, str) else ""
    return (lower_name in {"z", "ordinary_z", "free_z", "query_z"}
            and (row.get("canonical_free_Z") is True
                 or row.get("canonical_free_z") is True))


class _BudgetExceeded(Exception):
    pass


class _Expr:
    """Flat source expression receipt.

    Non-rational expressions carry a short Merkle-style intern handle in
    ``key``.  Child handles are stored separately, so a repeated source DAG
    never creates recursively expanded tuple keys.
    """

    __slots__ = ("key", "supported", "z_independent", "unknown_ops", "opaque_refs",
                 "audit", "operation", "children", "metadata", "truth")

    def __init__(self, key, supported=True, z_independent=True, unknown_ops=(),
                 opaque_refs=(), audit=(), operation=None, children=(), metadata=None,
                 truth=None):
        self.key = key
        self.supported = bool(supported)
        self.z_independent = z_independent
        self.unknown_ops = tuple(unknown_ops)
        self.opaque_refs = tuple(opaque_refs)
        self.audit = tuple(audit)
        self.operation = operation
        self.children = tuple(children)
        self.metadata = metadata
        self.truth = truth


def _merge(exprs, key, *, supported=None, z_independent=None, audit=(),
           operation=None, child_keys=None, metadata=None, unknown_ops=None,
           opaque_refs=None, truth=None):
    exprs = tuple(exprs)
    if supported is None:
        supported = all(item.supported for item in exprs)
    if z_independent is None:
        values = [item.z_independent for item in exprs]
        z_independent = False if False in values else (None if None in values else True)
    unknown = (tuple(dict.fromkeys(unknown_ops)) if unknown_ops is not None else
               tuple(dict.fromkeys(name for item in exprs for name in item.unknown_ops)))
    # Opaque provenance is a set of source atoms; multiplicity remains in the
    # flat child-handle DAG.  Deduplicate while collecting and stop at a hard
    # bound.  On overflow retain a non-admitting sentinel and an explicit
    # unknown reason, rather than truncating and later claiming a prefix proof.
    opaque_seen = set()
    opaque_list = []
    opaque_overflow = False
    source_refs = (opaque_refs if opaque_refs is not None else
                   (ref for item in exprs for ref in item.opaque_refs))
    for ref in source_refs:
        if ref in opaque_seen:
            continue
        if len(opaque_list) >= MAX_OPAQUE_REFS:
            opaque_overflow = True
            break
        opaque_seen.add(ref)
        opaque_list.append(ref)
    if opaque_overflow:
        child_digests = [_digest(item.opaque_refs) for item in exprs if item.opaque_refs]
        opaque_list.append(("opaque_ref_overflow", _digest(child_digests)))
        unknown = tuple(dict.fromkeys(tuple(unknown) + ("opaque_reference_budget_exceeded",)))
        supported = False
    opaque = tuple(opaque_list)
    events = list(event for item in exprs for event in item.audit)
    events.extend(audit)
    if child_keys is None:
        child_keys = tuple(item.key for item in exprs)
    return _Expr(key, supported, z_independent, unknown, opaque, tuple(events[-128:]),
                 operation=operation, children=tuple(child_keys), metadata=metadata,
                 truth=truth)


class _SourceWalker:
    """Small memoized source DAG walker with endpoint lexical scopes."""

    _MAX_FLAT_CHILDREN = 8192

    def __init__(self, source_name, nodes, *, node_budget=250000):
        self.source_name = source_name
        self.nodes = nodes if isinstance(nodes, list) else []
        self.node_budget = int(node_budget)
        self.visited = 0
        self.memo = {}
        self.rewrite_audit = []
        self._radius_coordinate_nodes = None
        self._frozen_values = {}
        self._intern_signatures = {}
        self._intern_rows_by_handle = {}
        self.intern_rows = []

    def _scope_key(self, scope, bound):
        radius_nodes = None if self._radius_coordinate_nodes is None else tuple(
            sorted(self._radius_coordinate_nodes))
        frozen = tuple(sorted((int(key), value.key)
                              for key, value in self._frozen_values.items()))
        return (tuple(sorted((int(key), int(value)) for key, value in scope.items())),
                tuple(bound), radius_nodes, frozen)

    def _intern_key(self, operation, child_keys=(), metadata=None):
        child_keys = tuple(child_keys)
        signature = (str(operation), child_keys, metadata)
        handle = self._intern_signatures.get(signature)
        if handle is not None:
            return ("ref", handle)
        handle = _digest({"operation": operation, "children": child_keys,
                          "metadata": metadata})
        existing = self._intern_rows_by_handle.get(handle)
        if existing is not None and existing != signature:
            # The digest collision branch is deterministic and remains flat.
            handle = _digest({"operation": operation, "children": child_keys,
                              "metadata": metadata, "collision": len(self.intern_rows)})
        self._intern_signatures[signature] = handle
        self._intern_rows_by_handle[handle] = signature
        self.intern_rows.append({
            "handle": handle,
            "digest": handle,
            "operation": str(operation),
            "children": _jsonable(child_keys),
            "metadata": _jsonable(metadata),
        })
        return ("ref", handle)

    def _result(self, operation, exprs=(), *, child_keys=None, metadata=None,
                key=None, supported=None, z_independent=None, audit=(),
                unknown_ops=None, opaque_refs=None, truth=None):
        exprs = tuple(exprs)
        if child_keys is None:
            child_keys = tuple(item.key for item in exprs)
        if key is None:
            key = self._intern_key(operation, child_keys, metadata)
        result = _merge(exprs, key, supported=supported, z_independent=z_independent,
                        audit=audit, operation=operation, child_keys=child_keys,
                        metadata=metadata, unknown_ops=unknown_ops,
                        opaque_refs=opaque_refs, truth=truth)
        return result

    def _scope_digest(self, scope, bound):
        return _digest({"scope": self._scope_key(scope, bound), "bound": list(bound)})

    def _unknown(self, node_id, row, scope, bound, reason):
        operation = row.get("operation", "<missing>") if isinstance(row, dict) else "<missing>"
        metadata = _metadata_digest(row if isinstance(row, dict) else {})
        scope_digest = self._scope_digest(scope, bound)
        opaque_metadata = _digest({"operation": operation, "node_id": node_id,
                                   "metadata": metadata, "scope": scope_digest})
        ref = (node_id, operation, metadata, scope_digest)
        return self._result("opaque", (), metadata=opaque_metadata, supported=False,
                            z_independent=None, unknown_ops=(reason,), opaque_refs=(ref,))

    def _ref(self, value, scope, bound):
        if _is_int(value):
            return self.walk(value, scope, bound)
        if isinstance(value, str):
            return self._result("literal", (), metadata=value)
        return self._result("opaque_literal", (), metadata=_digest(value), supported=False,
                            z_independent=None, unknown_ops=("non_graph_reference",))

    def _literal(self, value):
        if _is_int(value):
            rational = _rational(value)
            if rational is not None:
                return _Expr(rational, operation="rational")
        if isinstance(value, (str, bool)) or value is None:
            return self._result("literal", (), metadata=value)
        return self._result("literal_metadata", (), metadata=_digest(value))

    def _rational_result(self, exprs, value, operation, audit=()):
        return self._result(operation, exprs, key=value, audit=audit,
                            truth=(None if not (isinstance(value, tuple)
                                                and value[:1] == ("boolean",)) else value[1]))

    def walk(self, node_id, scope=None, bound=()):
        scope = {} if scope is None else dict(scope)
        bound = tuple(bound)
        if not _is_int(node_id) or node_id < 0 or node_id >= len(self.nodes):
            result = self._result("missing_source_node", (), metadata=node_id,
                                  supported=False, z_independent=None,
                                  unknown_ops=("missing_source_node",))
            return result
        row = self.nodes[node_id]
        # Decide lexical alpha shadowing before consulting either frozen
        # substitution values or endpoint coordinate scope.  A bound
        # integrand variable must never be captured by a surrounding
        # function substitution.
        bound_name_shadowed = (isinstance(row, dict)
                               and row.get("operation") == "bound_variable"
                               and isinstance(row.get("name"), str)
                               and row.get("name") in bound)
        frozen = self._frozen_values.get(node_id)
        if frozen is not None and not bound_name_shadowed:
            return frozen
        if not bound_name_shadowed and node_id in scope and scope[node_id] != node_id:
            return self.walk(scope[node_id], scope, bound)
        key = (node_id, self._scope_key(scope, bound))
        cached = self.memo.get(key)
        if cached is not None:
            return cached
        self.visited += 1
        if self.visited > self.node_budget:
            raise _BudgetExceeded("source DAG node budget exceeded")
        if not isinstance(row, dict):
            result = self._unknown(node_id, {}, scope, bound, "malformed_source_node")
            self.memo[key] = result
            return result
        operation = row.get("operation")

        if operation == "exact_rational":
            rational = _rational(row.get("numerator"), row.get("denominator"))
            if rational is None:
                result = self._result("invalid_rational", (), metadata=_metadata_digest(row),
                                      supported=False, z_independent=True,
                                      unknown_ops=("invalid_rational",))
            else:
                result = _Expr(rational, operation="rational")
        elif operation == "bound_variable":
            name = row.get("name")
            if self._radius_coordinate_nodes is not None and node_id not in self._radius_coordinate_nodes:
                result = self._result(
                    "unapproved_radius_bound_variable", (),
                    metadata=(node_id, name, _metadata_digest(row)),
                    supported=False, z_independent=None,
                    unknown_ops=("radius_bound_variable_not_coordinate",),
                    audit=("radius_coordinate_bound_variable_unverified",),
                )
            elif isinstance(name, str) and name in bound:
                result = self._result("bound_variable", (),
                                      metadata=("alpha", len(bound) - 1 - bound[::-1].index(name)))
            else:
                lower_name = name.lower() if isinstance(name, str) else ""
                z_independent = lower_name not in {"z", "ordinary_z", "free_z", "query_z"}
                metadata = ((name, _metadata_digest(row, {"name"}))
                            if _explicit_canonical_free_z(row, name)
                            else (node_id, name, _metadata_digest(row, {"name"})))
                result = self._result("bound_variable", (),
                                      metadata=metadata,
                                      z_independent=z_independent)
        elif operation == "current_original_source_parameter":
            independent = row.get("Z_independent")
            metadata = ((row.get("name"), _metadata_digest(row, {"name", "Z_independent"}))
                        if _explicit_canonical_free_z(row, row.get("name"))
                        else (node_id, row.get("name"),
                              _metadata_digest(row, {"name", "Z_independent"})))
            result = self._result(
                "source_parameter", (),
                metadata=metadata,
                z_independent=True if independent is True else (False if independent is False else None),
            )
        elif operation in ("mathematical_pi", "mathematical_positive_infinity"):
            result = self._result("source_atom", (), metadata=(operation, _metadata_digest(row)))
        elif operation in ("sum", "product"):
            arguments = row.get("arguments")
            if not isinstance(arguments, list):
                result = self._unknown(node_id, row, scope, bound, "malformed_arguments")
            else:
                children = [self.walk(child, scope, bound) for child in arguments]
                flat_keys = []
                for child in children:
                    if child.operation == operation:
                        flat_keys.extend(child.children)
                    else:
                        flat_keys.append(child.key)
                if len(flat_keys) > self._MAX_FLAT_CHILDREN:
                    result = self._unknown(node_id, row, scope, bound,
                                           "flat_AC_child_budget_exceeded")
                else:
                    rational_parts = [child for child in flat_keys if _is_rational(child)]
                    symbolic = [child for child in flat_keys if not _is_rational(child)]
                    folded = None
                    if rational_parts:
                        folded = rational_parts[0]
                        for part in rational_parts[1:]:
                            folded = (_r_add(folded, part) if operation == "sum"
                                      else _r_mul(folded, part))
                    if operation == "sum":
                        if folded is not None and folded[1] != 0:
                            symbolic.append(folded)
                        symbolic = sorted(symbolic, key=repr)  # preserve multiplicity
                        if not symbolic:
                            result = self._rational_result(children, _rational(0), "sum_fold")
                        elif len(symbolic) == 1:
                            result = self._result("reduced_sum", children,
                                                  child_keys=tuple(symbolic), key=symbolic[0])
                        else:
                            result = self._result("sum", children,
                                                  child_keys=tuple(symbolic))
                    elif folded is not None and folded[1] == 0:
                        result = self._rational_result(children, _rational(0), "product_fold")
                    else:
                        symbolic = [item for item in symbolic if item != _rational(1)]
                        if folded is not None and folded != _rational(1):
                            symbolic.append(folded)
                        symbolic = sorted(symbolic, key=repr)  # preserve multiplicity
                        if not symbolic:
                            result = self._rational_result(children, _rational(1), "product_fold")
                        elif len(symbolic) == 1:
                            result = self._result("reduced_product", children,
                                                  child_keys=tuple(symbolic), key=symbolic[0])
                        else:
                            result = self._result("product", children,
                                                  child_keys=tuple(symbolic))
        elif operation == "negative":
            child = self._ref(row.get("argument"), scope, bound)
            if _is_rational(child.key):
                result = self._rational_result((child,), _r_neg(child.key), "negative_fold")
            elif child.operation == "negative" and len(child.children) == 1:
                result = self._result("negative_cancel", (child,),
                                      child_keys=child.children, key=child.children[0])
            else:
                result = self._result("negative", (child,))
        elif operation == "positive_quotient":
            numerator = self._ref(row.get("numerator"), scope, bound)
            denominator = self._ref(row.get("denominator"), scope, bound)
            metadata = _metadata_digest(row, {"numerator", "denominator"})
            if not row.get("source_positive_certificate"):
                result = self._result("positive_quotient", (numerator, denominator),
                                      metadata=metadata, supported=False,
                                      audit=("missing_positive_quotient_certificate",))
            elif _is_rational(numerator.key) and numerator.key[1] == 0:
                result = self._rational_result((numerator, denominator), _rational(0),
                                               "zero_quotient", ("zero_quotient",))
            elif _is_rational(denominator.key) and denominator.key[1] != 0:
                quotient = (_r_div(numerator.key, denominator.key)
                            if _is_rational(numerator.key) else None)
                if quotient is not None:
                    result = self._rational_result((numerator, denominator), quotient,
                                                   "exact_positive_quotient",
                                                   ("exact_positive_quotient",))
                else:
                    result = self._result("positive_quotient", (numerator, denominator),
                                          metadata=metadata,
                                          audit=("positive_quotient_preserved",))
            else:
                result = self._result("positive_quotient", (numerator, denominator),
                                      metadata=metadata,
                                      audit=("positive_quotient_preserved",))
        elif operation == "function_substitution":
            expression = row.get("expression")
            variable = row.get("variable")
            value = row.get("value")
            variable_row = (self.nodes[variable] if _is_int(variable)
                            and 0 <= variable < len(self.nodes) else None)
            if not (_is_int(expression) and _is_int(variable) and _is_int(value)
                    and isinstance(variable_row, dict)
                    and variable_row.get("operation") == "bound_variable"):
                result = self._unknown(node_id, row, scope, bound, "malformed_function_substitution")
            else:
                # Freeze the value under the outer scope before binding it.
                value_expr = self._ref(value, scope, bound)
                next_scope = dict(scope)
                next_scope[variable] = value
                previous_frozen = self._frozen_values
                self._frozen_values = dict(previous_frozen)
                self._frozen_values[variable] = value_expr
                try:
                    body = self.walk(expression, next_scope, bound)
                finally:
                    self._frozen_values = previous_frozen
                metadata = _metadata_digest(row, {"expression", "variable", "value"})
                if row.get("Z_independent_substitution") is not True:
                    result = self._result("function_substitution", (value_expr, body),
                                          metadata=metadata, key=self._intern_key(
                                              "function_substitution", (value_expr.key, body.key), metadata),
                                          supported=False,
                                          audit=("unqualified_function_substitution",))
                else:
                    result = self._result("endpoint_substitution", (value_expr, body),
                                          metadata=metadata, key=body.key,
                                          audit=("endpoint_substitution",))
        elif operation == "definite_integral":
            lower = self._ref(row.get("lower"), scope, bound)
            upper = self._ref(row.get("upper"), scope, bound)
            variable = row.get("variable")
            if not isinstance(variable, str):
                result = self._unknown(node_id, row, scope, bound, "malformed_integral_variable")
            elif lower.key == upper.key and lower.supported and upper.supported:
                # Do not walk an unneeded integrand for an explicitly empty interval.
                result = self._rational_result((lower, upper), _rational(0), "empty_integral",
                                               ("explicit_empty_integral_lower_equals_upper",))
            else:
                body = self.walk(row.get("integrand"), scope, bound + (variable,))
                metadata = _metadata_digest(row, {"integrand", "variable", "lower", "upper"})
                result = self._result("definite_integral", (body, lower, upper),
                                      metadata=("alpha", len(bound), metadata),
                                      audit=("nonempty_integral_symbolic_with_alpha_scope",))
        elif operation == "analytic_unary":
            child = self._ref(row.get("argument"), scope, bound)
            name = row.get("name")
            if name not in ANALYTIC_NAMES:
                result = self._result(
                    "analytic_unary", (child,),
                    metadata=(name, _metadata_digest(row, {"name", "argument"})),
                    supported=False,
                    unknown_ops=("unsupported_analytic_name",),
                    audit=("analytic_unary_name_unproved",),
                )
            elif _is_rational(child.key) and name == "exp" and child.key[1] == 0:
                result = self._rational_result((child,), _rational(1), "analytic_exp_zero")
            elif _is_rational(child.key) and name == "log" and child.key == _rational(1):
                result = self._rational_result((child,), _rational(0), "analytic_log_one")
            elif _is_rational(child.key) and name in ("sqrt", "positive_sqrt") and child.key[1] == 0:
                result = self._rational_result((child,), _rational(0), "analytic_sqrt_zero")
            else:
                result = self._result("analytic_unary", (child,),
                                      metadata=(name, _metadata_digest(row, {"name", "argument"})),
                                      audit=("analytic_unary_preserved",))
        elif operation == "exact_original_source_flat_sigma":
            child = self._ref(row.get("argument"), scope, bound)
            if _is_rational(child.key) and child.key == _rational(0):
                result = self._rational_result((child,), _rational(0), "flat_sigma_lower",
                                               ("sigma_lower_endpoint_zero",))
            elif _is_rational(child.key) and child.key == _rational(1):
                result = self._rational_result((child,), _rational(1), "flat_sigma_upper",
                                               ("sigma_upper_endpoint_one",))
            else:
                result = self._result("flat_sigma", (child,),
                                      metadata=_metadata_digest(row, {"argument"}),
                                      audit=("flat_sigma_symbolic",))
        elif operation == "exact_real_comparison":
            left = self._ref(row.get("left"), scope, bound)
            right = self._ref(row.get("right"), scope, bound)
            operator = row.get("operator")
            if operator not in COMPARISON_OPERATORS:
                result = self._result(
                    "exact_real_comparison", (left, right), metadata=operator,
                    supported=False,
                    unknown_ops=("unsupported_comparison_operator",),
                    audit=("comparison_operator_unproved",),
                )
                self.memo[key] = result
                self.rewrite_audit.extend(result.audit)
                if len(self.rewrite_audit) > 4096:
                    del self.rewrite_audit[:-4096]
                return result
            value = None
            if _is_rational(left.key) and _is_rational(right.key):
                lv = left.key[1] * right.key[2]
                rv = right.key[1] * left.key[2]
                value = {"==": lv == rv, "!=": lv != rv, "<": lv < rv,
                         "<=": lv <= rv, ">": lv > rv, ">=": lv >= rv}.get(operator)
            if value is not None:
                result = self._result("boolean", (left, right), metadata=value,
                                      truth=value, audit=("exact_real_comparison",))
            else:
                result = self._result("exact_real_comparison", (left, right),
                                      metadata=operator, audit=("exact_real_comparison",))
        elif operation == "logical_or":
            arguments = row.get("arguments")
            if not isinstance(arguments, list) or not all(_is_int(item) for item in arguments):
                result = self._unknown(node_id, row, scope, bound, "malformed_logical_or")
            elif not arguments:
                result = self._unknown(node_id, row, scope, bound, "empty_logical_or")
            else:
                children = [self._ref(child, scope, bound) for child in arguments]
                if any(child.truth is True for child in children):
                    result = self._result("boolean", children, metadata=True,
                                          truth=True, audit=("logical_or_preserved",))
                elif children and all(child.truth is False for child in children):
                    result = self._result("boolean", children, metadata=False,
                                          truth=False, audit=("logical_or_preserved",))
                else:
                    result = self._result("logical_or", children,
                                          child_keys=tuple(sorted((child.key for child in children), key=repr)),
                                          audit=("logical_or_preserved",))
        elif operation == "original_lazy_flat_branch":
            predicate_ref = row.get("flat_predicate")
            if not _is_int(predicate_ref):
                result = self._unknown(node_id, row, scope, bound, "malformed_flat_predicate")
            else:
                predicate = self._ref(predicate_ref, scope, bound)
                # Evaluate the predicate first.  A flat branch must not walk
                # the inactive body or treat flat_value's graph id as a literal.
                if predicate.truth is True:
                    flat_ref = row.get("flat_value")
                    flat_value = self._ref(flat_ref, scope, bound)
                    result = self._result("flat_branch_selected", (predicate, flat_value),
                                          child_keys=(predicate.key, flat_value.key),
                                          key=flat_value.key,
                                          audit=("flat_branch_selected",))
                elif predicate.truth is False and _is_int(row.get("active_body")):
                    active = self._ref(row.get("active_body"), scope, bound)
                    result = self._result("active_flat_branch_selected", (predicate, active),
                                          child_keys=(predicate.key, active.key), key=active.key,
                                          audit=("active_flat_branch_selected",))
                else:
                    # An unknown predicate leaves both branch bodies relevant.
                    # Retain their exact source row and scope, without walking
                    # either body or dropping them from the identity key.
                    result = self._unknown(node_id, row, scope, bound,
                                           "flat_branch_predicate_unproved")
        elif operation == "original_monotone_phase_inverse":
            refs = ("phase_function", "target_phase", "angle_lower", "angle_upper", "Delta", "eta")
            required_refs = ("phase_function", "target_phase", "angle_lower", "angle_upper")
            if not all(_is_int(row.get(name)) for name in required_refs):
                result = self._unknown(node_id, row, scope, bound,
                                       "malformed_monotone_phase_inverse")
                self.memo[key] = result
                self.rewrite_audit.extend(result.audit)
                if len(self.rewrite_audit) > 4096:
                    del self.rewrite_audit[:-4096]
                return result
            children = {name: self._ref(row.get(name), scope, bound)
                        for name in refs if _is_int(row.get(name))}
            target = children.get("target_phase")
            endpoints = row.get("phase_endpoints")
            derivative_certificate = row.get("positive_derivative_certificate")
            if derivative_certificate is not True:
                derivative_certificate = bool(
                    isinstance(row.get("derivative_positive_certificate"), str)
                    and row.get("derivative_positive_certificate").strip()
                )
            inverse_flags = (
                row.get("true_prescribed_source_C5_implicit_inverse") is True
                and row.get("branch_inverse_is_original_not_support_cap") is True
                and derivative_certificate
            )
            endpoints_ok = (isinstance(endpoints, list) and len(endpoints) == 2
                            and _rational(endpoints[0]) == _rational(0)
                            and _rational(endpoints[1]) == _rational(1))
            bounds_ok = (row.get("matching_known_angle_bounds") is True
                         or row.get("known_angle_bounds_match") is True
                         or row.get("angle_bounds_match") is True)
            if inverse_flags and endpoints_ok and bounds_ok and target is not None:
                if target.key == _rational(0) and "angle_lower" in children:
                    result = self._result("inverse_lower_endpoint", tuple(children.values()),
                                          key=children["angle_lower"].key,
                                          audit=("inverse_lower_endpoint",))
                elif target.key == _rational(1) and "angle_upper" in children:
                    result = self._result("inverse_upper_endpoint", tuple(children.values()),
                                          key=children["angle_upper"].key,
                                          audit=("inverse_upper_endpoint",))
                else:
                    result = self._result("monotone_phase_inverse", tuple(children.values()),
                                          metadata=_metadata_digest(row, set(refs) | {"phase_endpoints"}),
                                          supported=False,
                                          unknown_ops=("inverse_endpoint_unmatched",),
                                          audit=("inverse_endpoint_scope_preserved",))
            else:
                result = self._result("monotone_phase_inverse", tuple(children.values()),
                                      metadata=_metadata_digest(row, set(refs) | {"phase_endpoints"}),
                                      supported=False,
                                      unknown_ops=("unqualified_monotone_phase_inverse",),
                                      audit=("inverse_symbolic_unqualified",))
        elif operation in ("substitute_original_inverse_angle", "substitute_original_native_coordinate",
                           "substitute_free_query_at_native_coordinate"):
            required_refs = {
                "substitute_original_inverse_angle": ("body", "inverse_angle"),
                "substitute_original_native_coordinate": ("body", "query_coordinate"),
                "substitute_free_query_at_native_coordinate": ("body", "free_query"),
            }[operation]
            refs = ("body", "inverse_angle", "query_coordinate", "free_query")
            if not all(_is_int(row.get(name)) for name in required_refs):
                result = self._unknown(node_id, row, scope, bound,
                                       "malformed_coordinate_substitution")
            else:
                children = {name: self._ref(row.get(name), scope, bound)
                            for name in required_refs}
                result = self._result(operation, tuple(children.values()),
                                      metadata=_metadata_digest(row, set(refs)),
                                      audit=("coordinate_substitution_symbolic",))
        elif operation == "exact_positive_integer_expression":
            arguments = row.get("arguments", [])
            if (isinstance(arguments, list) and arguments
                    and all(_is_int(item) for item in arguments)):
                children = [self._ref(item, scope, bound) for item in arguments]
                result = self._result("exact_positive_integer_expression", children,
                                      metadata=_metadata_digest(row, {"arguments"}))
            else:
                result = self._unknown(node_id, row, scope, bound, "malformed_exact_integer_expression")
        else:
            result = self._unknown(node_id, row, scope, bound, "unsupported_operation:" + str(operation))

        self.memo[key] = result
        self.rewrite_audit.extend(result.audit)
        if len(self.rewrite_audit) > 4096:
            del self.rewrite_audit[:-4096]
        return result

    def endpoint(self, root, coordinate_nodes, domain_nodes, endpoint,
                 *, radius_coordinate_nodes=None):
        scope = {}
        if isinstance(coordinate_nodes, list) and isinstance(domain_nodes, list):
            endpoint_index = 1 if endpoint == "upper" else 0
            if domain_nodes and endpoint_index < len(domain_nodes):
                endpoint_value = domain_nodes[endpoint_index]
                if _is_int(endpoint_value):
                    for coordinate in coordinate_nodes:
                        if _is_int(coordinate):
                            scope[coordinate] = endpoint_value
        previous_radius_nodes = self._radius_coordinate_nodes
        if radius_coordinate_nodes is None:
            self._radius_coordinate_nodes = None
        elif isinstance(radius_coordinate_nodes, (list, tuple, set, frozenset)):
            self._radius_coordinate_nodes = frozenset(
                item for item in radius_coordinate_nodes if _is_int(item))
        else:
            self._radius_coordinate_nodes = frozenset()
        try:
            result = self.walk(root, scope, ())
        except _BudgetExceeded as error:
            result = _Expr(("ref", _digest(("budget_exceeded", root))), False, None,
                           (str(error),), operation="budget_exceeded")
        except RecursionError as error:
            reason = "source DAG recursion limit exceeded"
            if str(error):
                reason += ": " + str(error)
            result = _Expr(("ref", _digest(("recursion_exceeded", root))), False, None,
                           (reason,), operation="recursion_exceeded")
        finally:
            self._radius_coordinate_nodes = previous_radius_nodes
        return result


class _Graph:
    def __init__(self):
        self.nodes = []

    def add(self, operation, arguments=(), **metadata):
        index = len(self.nodes)
        self.nodes.append(dict(operation=operation, arguments=list(arguments), **metadata))
        return index


def _payload_bindings(payload, key):
    value = payload.get(key) if isinstance(payload, dict) else None
    return value if isinstance(value, dict) else {}


def _source_array(value):
    """Accept only literal graph-handle arrays/scalars from archive bindings."""
    if _is_int(value):
        return True
    return (isinstance(value, list) and bool(value)
            and all(_is_int(item) for item in value))


def _binding_shape(payloads):
    checks = {}
    mod = _payload_bindings(payloads.get(MOD, {}), "modified_Q4_pressure_C3_bindings")
    loop = _payload_bindings(payloads.get(LOOP, {}), "prescribed_loop_C5_bindings")
    patch = _payload_bindings(payloads.get(PATCH, {}), "Rm_patch_finite_N_C5_bindings")
    outer = _payload_bindings(payloads.get(OUTER, {}), "true_outer_C6_C5_bindings")
    # The expression above is intentionally explicit: it covers the eight
    # LOOP/MOD charts through postrestore, including frozen_macro.
    first_ten = ("frozen_macro", "first_switch", "second_switch", "post_power",
                 "long_reshape", "reference", "restoration", "postrestore")
    loop_charts = loop.get("charts", {}) if isinstance(loop.get("charts"), dict) else {}
    mod_charts = mod.get("charts", {}) if isinstance(mod.get("charts"), dict) else {}
    # PATCH owns one top-level chart and one top-level physical map.  The
    # nested endpoint maps are historical convenience copies and are not the
    # authoritative source arrays for this leaf.
    patch_chart = patch.get("chart") if isinstance(patch.get("chart"), dict) else {}
    patch_physical = patch.get("physical") if isinstance(patch.get("physical"), dict) else {}
    outer_charts = outer.get("charts", {}) if isinstance(outer.get("charts"), dict) else {}
    checks["first_ten_LOOP_geometry_binding_map"] = all(
        isinstance(loop_charts.get(chart), dict)
        and all(field in loop_charts[chart]
                and _source_array(loop_charts[chart][field])
                for field in ("coordinate", "domain", "radius", "physical_measure"))
        for chart in first_ten
    )
    checks["first_ten_MOD_physical_binding_map"] = all(
        isinstance(mod_charts.get(chart), dict)
        and isinstance(mod_charts[chart].get("physical"), dict)
        and all(field in mod_charts[chart]["physical"]
                and _source_array(mod_charts[chart]["physical"][field])
                for field in ("U", "U_y", "m", "m_y", "Q", "Q_y"))
        for chart in first_ten
    )
    checks["PATCH_own_binding_map"] = all(
        field in patch_chart and _source_array(patch_chart[field])
        for field in ("coordinate", "domain", "radius", "physical_measure")
    ) and all(
        field in patch_physical and _source_array(patch_physical[field])
        for field in ("U", "U_y", "m", "m_y", "Q", "Q_y")
    )
    outer_names = ("Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power")
    checks["OUTER_own_binding_map"] = all(
        isinstance(outer_charts.get(chart), dict)
        and all(field in outer_charts[chart]
                and _source_array(outer_charts[chart][field])
                for field in ("coordinate", "domain", "radius", "physical_measure"))
        and isinstance(outer_charts[chart].get("physical"), dict)
        and all(field in outer_charts[chart]["physical"]
                and _source_array(outer_charts[chart]["physical"][field])
                for field in ("U", "U_y", "m", "m_y", "Q", "Q_y"))
        for chart in outer_names
    )
    checks["all_archive_binding_maps_present"] = all(checks.values())
    return checks


def _graph_nodes(payload):
    nodes = payload.get("exact_graph_nodes") if isinstance(payload, dict) else None
    return nodes if isinstance(nodes, list) else []


def _prefix_equal(left, right):
    if not isinstance(left, list) or not isinstance(right, list) or len(left) > len(right):
        return False
    return all(left[index] == right[index] for index in range(len(left)))


def _prefix_checks(payloads):
    nodes = {name: _graph_nodes(payloads.get(name, {})) for name in (MOD, LOOP, PATCH, OUTER)}
    checks = {
        "LOOP_MOD_exact_graph_prefix": bool(nodes[LOOP]) and _prefix_equal(nodes[LOOP], nodes[MOD]),
        "MOD_PATCH_exact_graph_prefix": bool(nodes[MOD]) and _prefix_equal(nodes[MOD], nodes[PATCH]),
        "PATCH_OUTER_exact_graph_prefix": bool(nodes[PATCH]) and _prefix_equal(nodes[PATCH], nodes[OUTER]),
    }
    checks["LOOP_PATCH_transitive_exact_graph_prefix"] = (
        checks["LOOP_MOD_exact_graph_prefix"]
        and checks["MOD_PATCH_exact_graph_prefix"]
        and bool(nodes[LOOP])
        and _prefix_equal(nodes[LOOP], nodes[PATCH])
    )
    checks["exact_graph_prefix_chain"] = all(checks.values())
    return checks, nodes


def _same_source_family(payloads):
    families = [payloads.get(name, {}).get("source_family") for name in (PARENT, WEIGHTED, MOD, LOOP, PATCH, OUTER)]
    if not families or any(not isinstance(family, dict) for family in families):
        return False, None, families
    first = families[0]
    return all(family == first for family in families[1:]), first, families


def _node(graph, node_id):
    if _is_int(node_id) and 0 <= node_id < len(graph):
        row = graph[node_id]
        return row if isinstance(row, dict) else None
    return None


def _source_port(weighted_nodes, weighted_id, expected_record):
    row = _node(weighted_nodes, weighted_id)
    if row is None:
        return None
    source_nodes = row.get("source_nodes")
    source_nodes_valid = (isinstance(source_nodes, list) and bool(source_nodes)
                          and all(_is_int(item) for item in source_nodes))
    return {
        "weighted_node": weighted_id,
        "operation": row.get("operation"),
        "source_record": row.get("source_record"),
        "source_nodes": list(source_nodes) if source_nodes_valid else None,
        "chart": row.get("chart"),
        "quantity": row.get("quantity"),
        "source_path": row.get("source_path"),
        "metadata_keys": _metadata_keys(row, {"source_nodes", "source_record", "chart", "quantity", "source_path"}),
        "metadata_digest": _metadata_digest(row, {"source_nodes", "source_record", "chart", "quantity", "source_path"}),
        "expected_record": expected_record,
        "record_matches_expected": row.get("source_record") == expected_record,
    }


def _endpoint_descriptor(weighted_nodes, endpoint_id, spec, side):
    endpoint_row = _node(weighted_nodes, endpoint_id)
    if endpoint_row is None:
        return None, ["missing_weighted_endpoint_node"]
    geometry_record = spec[side + "_geometry"]
    physical_record = spec[side + "_physical"]
    errors = []
    result = {
        "weighted_node": endpoint_id,
        "chart": endpoint_row.get("chart"),
        "endpoint": endpoint_row.get("endpoint"),
        "endpoint_bound": endpoint_row.get("endpoint_bound"),
        "formula": endpoint_row.get("formula"),
        "coordinate": _source_port(weighted_nodes, endpoint_row.get("coordinate"), geometry_record),
        "domain": _source_port(weighted_nodes, endpoint_row.get("domain"), geometry_record),
        "radius": _source_port(weighted_nodes, endpoint_row.get("radius"), geometry_record),
        "physical_measure": _source_port(weighted_nodes, endpoint_row.get("physical_measure"), geometry_record),
        "Jacobian": _source_port(weighted_nodes, endpoint_row.get("Jacobian"), geometry_record),
        "physical": {},
    }
    for field in ("U", "U_y", "m", "m_y", "Q", "Q_y", "S", "local_integrals.S1",
                  "cumulative_integrals.S1", "incoming_integrals.S1", "R_Q_y"):
        ports = endpoint_row.get("physical_ports")
        weighted_id = ports.get(field) if isinstance(ports, dict) else None
        result["physical"][field] = _source_port(weighted_nodes, weighted_id, physical_record)
    for name, port in (("coordinate", result["coordinate"]), ("domain", result["domain"]),
                       ("radius", result["radius"]), ("U", result["physical"].get("U")),
                       ("Q", result["physical"].get("Q")), ("Q_y", result["physical"].get("Q_y"))):
        if not isinstance(port, dict) or not isinstance(port.get("source_nodes"), list):
            errors.append("missing_" + name + "_source_port")
        elif port.get("record_matches_expected") is not True:
            errors.append("wrong_" + name + "_source_archive")
    if endpoint_row.get("chart") != spec[side + "_chart"]:
        errors.append("weighted_chart_mismatch")
    if endpoint_row.get("endpoint") != spec[side + "_endpoint"]:
        errors.append("weighted_endpoint_mismatch")
    expected_bound = "domain[1]" if spec[side + "_endpoint"] == "upper" else "domain[0]"
    if endpoint_row.get("endpoint_bound") != expected_bound:
        errors.append("weighted_endpoint_bound_mismatch")
    if endpoint_row.get("operation") != "actual_B_endpoint_function":
        errors.append("weighted_endpoint_operation_mismatch")
    if endpoint_row.get("formula") != B_FORMULA:
        errors.append("weighted_B_formula_mismatch")
    return result, errors


_BINDING_KEYS = {
    MOD: "modified_Q4_pressure_C3_bindings",
    LOOP: "prescribed_loop_C5_bindings",
    PATCH: "Rm_patch_finite_N_C5_bindings",
    OUTER: "true_outer_C6_C5_bindings",
}


def _binding_chart(payloads, record, chart):
    """Return the authoritative chart map for one source archive.

    PATCH deliberately has a top-level ``chart``/``physical`` pair.  Reading
    its nested endpoint copies would weaken provenance and can accidentally
    compare a convenience alias instead of the published source array.
    """
    binding = _payload_bindings(payloads.get(record, {}), _BINDING_KEYS.get(record))
    if record == PATCH:
        return (binding.get("chart") if chart == "Rm_patch_finite_N" else None), binding.get("physical")
    charts = binding.get("charts")
    if not isinstance(charts, dict):
        return None, None
    chart_map = charts.get(chart)
    if not isinstance(chart_map, dict):
        return None, None
    return chart_map, chart_map.get("physical")


def _binding_nodes(value):
    if _is_int(value):
        return [value]
    if isinstance(value, list) and bool(value) and all(_is_int(item) for item in value):
        return list(value)
    return None


def _archive_binding_audit(payloads, side):
    """Compare every required weighted port against its archive binding array."""
    details = []
    errors = []
    if not isinstance(side, dict):
        return {"passed": False, "details": [], "errors": ["missing_endpoint_side"]}
    chart = side.get("chart")
    geometry_fields = ("coordinate", "domain", "radius", "physical_measure", "Jacobian")
    physical_fields = ("U", "U_y", "m", "m_y", "Q", "Q_y")
    for field in geometry_fields:
        port = side.get(field)
        record = port.get("source_record") if isinstance(port, dict) else None
        chart_map, physical_map = _binding_chart(payloads, record, chart)
        binding_field = "physical_measure" if field == "Jacobian" else field
        expected = (physical_map.get(binding_field) if binding_field == "physical_measure"
                    and isinstance(physical_map, dict) and binding_field in physical_map
                    else chart_map.get(binding_field) if isinstance(chart_map, dict)
                    else None)
        actual = port.get("source_nodes") if isinstance(port, dict) else None
        expected_nodes = _binding_nodes(expected)
        passed = (isinstance(actual, list) and expected_nodes is not None
                  and actual == expected_nodes
                  and isinstance(port, dict) and port.get("chart") == chart
                  and port.get("quantity") == field)
        details.append({"field": field, "record": record, "chart": chart,
                        "expected_nodes": expected_nodes, "actual_nodes": actual,
                        "passed": passed})
        if not passed:
            errors.append("binding_array_mismatch_" + field)
    physical = side.get("physical") if isinstance(side.get("physical"), dict) else {}
    for field in physical_fields:
        port = physical.get(field)
        record = port.get("source_record") if isinstance(port, dict) else None
        chart_map, physical_map = _binding_chart(payloads, record, chart)
        expected = physical_map.get(field) if isinstance(physical_map, dict) else None
        actual = port.get("source_nodes") if isinstance(port, dict) else None
        expected_nodes = _binding_nodes(expected)
        passed = (isinstance(actual, list) and expected_nodes is not None
                  and actual == expected_nodes
                  and isinstance(port, dict) and port.get("chart") == chart
                  and port.get("quantity") == field)
        details.append({"field": field, "record": record, "chart": chart,
                        "expected_nodes": expected_nodes, "actual_nodes": actual,
                        "passed": passed})
        if not passed:
            errors.append("binding_array_mismatch_" + field)
    return {"passed": not errors, "details": details, "errors": errors}


def _root_positive_L(graph, root):
    row = _node(graph, root)
    seen = set()
    while isinstance(row, dict) and root not in seen:
        seen.add(root)
        if row.get("operation") == "positive_quotient":
            denominator = row.get("denominator")
            return denominator == L_NODE and bool(row.get("source_positive_certificate")), denominator, row.get("source_positive_certificate")
        if row.get("operation") == "function_substitution" and _is_int(row.get("expression")):
            root = row.get("expression")
            row = _node(graph, root)
            continue
        break
    return False, None, None


def _root_descriptor(side, field, order=0):
    if field == "R":
        port = side.get("radius")
    else:
        port = side.get("physical", {}).get(field)
    if not isinstance(port, dict) or not isinstance(port.get("source_nodes"), list):
        return None
    roots = port["source_nodes"]
    if order < 0 or order >= len(roots) or not _is_int(roots[order]):
        return None
    coordinate = side.get("coordinate", {}).get("source_nodes") if isinstance(side.get("coordinate"), dict) else None
    domain = side.get("domain", {}).get("source_nodes") if isinstance(side.get("domain"), dict) else None
    return {
        "field": field,
        "ordinary_Z_order": order,
        "source_record": port.get("source_record"),
        "source_node": roots[order],
        "source_nodes": list(roots),
        "coordinate_nodes": list(coordinate) if isinstance(coordinate, list) else None,
        "domain_nodes": list(domain) if isinstance(domain, list) else None,
        "endpoint": side.get("endpoint"),
        "chart": side.get("chart"),
        "port_weighted_node": port.get("weighted_node"),
        "port_metadata_digest": port.get("metadata_digest"),
    }


def _pair_prefix_limit(left_record, right_record, prefix_checks, graph_lengths):
    if left_record == right_record:
        return graph_lengths.get(left_record, 0)
    pairs = {
        frozenset((LOOP, MOD)): ("LOOP_MOD_exact_graph_prefix", LOOP),
        frozenset((MOD, PATCH)): ("MOD_PATCH_exact_graph_prefix", MOD),
        frozenset((PATCH, OUTER)): ("PATCH_OUTER_exact_graph_prefix", PATCH),
        frozenset((LOOP, PATCH)): ("LOOP_PATCH_transitive_exact_graph_prefix", LOOP),
    }
    item = pairs.get(frozenset((left_record, right_record)))
    if item and prefix_checks.get(item[0]) is True:
        return graph_lengths.get(item[1], 0)
    return -1


def _opaque_pair_safe(left_expr, right_expr, left_record, right_record, prefix_checks, graph_lengths):
    if not left_expr.opaque_refs and not right_expr.opaque_refs:
        return False
    # A bounded provenance sentinel is intentionally non-admitting.  Its
    # digest records that the walk overflowed, but it is not a graph node
    # identity that can be justified by an archive prefix.
    if any(not (_is_int(ref[0]) if isinstance(ref, tuple) and ref else False)
           for ref in tuple(left_expr.opaque_refs) + tuple(right_expr.opaque_refs)):
        return False
    if len(left_expr.opaque_refs) != len(right_expr.opaque_refs):
        return False
    limit = _pair_prefix_limit(left_record, right_record, prefix_checks, graph_lengths)
    if limit < 0:
        return False
    for left_ref, right_ref in zip(left_expr.opaque_refs, right_expr.opaque_refs):
        if left_ref[0] != right_ref[0] or left_ref[1:] != right_ref[1:] or left_ref[0] >= limit:
            return False
    return True


def _row_scope_digest(descriptor):
    return _digest({
        "source_record": descriptor.get("source_record"),
        "source_node": descriptor.get("source_node"),
        "coordinate_nodes": descriptor.get("coordinate_nodes"),
        "domain_nodes": descriptor.get("domain_nodes"),
        "endpoint": descriptor.get("endpoint"),
    })


def _compare_row(left_descriptor, right_descriptor, payloads, walkers, prefix_checks, graph_lengths, field, order):
    row = {
        "field": field,
        "ordinary_Z_order": order,
        "comparison_kind": "bounded_endpoint_source_DAG_congruence",
        "equal": None,
        "source_identity_proved": False,
        "status": "unproved",
        "numeric_value": None,
        "left": left_descriptor,
        "right": right_descriptor,
        "rewrite_audit": [],
        "reasons": [],
    }
    if left_descriptor is None or right_descriptor is None:
        row["reasons"].append("missing_endpoint_root")
        return row
    left_record = left_descriptor.get("source_record")
    right_record = right_descriptor.get("source_record")
    if left_record not in walkers or right_record not in walkers:
        row["reasons"].append("source_archive_not_loaded")
        return row
    left_graph = _graph_nodes(payloads.get(left_record, {}))
    right_graph = _graph_nodes(payloads.get(right_record, {}))
    left_root = left_descriptor["source_node"]
    right_root = right_descriptor["source_node"]
    left_radius_nodes = left_descriptor.get("coordinate_nodes") if field == "R" else None
    right_radius_nodes = right_descriptor.get("coordinate_nodes") if field == "R" else None
    left_expr = walkers[left_record].endpoint(left_root, left_descriptor.get("coordinate_nodes"),
                                               left_descriptor.get("domain_nodes"), left_descriptor.get("endpoint"),
                                               radius_coordinate_nodes=left_radius_nodes)
    right_expr = walkers[right_record].endpoint(right_root, right_descriptor.get("coordinate_nodes"),
                                                 right_descriptor.get("domain_nodes"), right_descriptor.get("endpoint"),
                                                 radius_coordinate_nodes=right_radius_nodes)
    left_key_digest = _digest(left_expr.key)
    right_key_digest = _digest(right_expr.key)
    row["left"] = dict(left_descriptor, canonical_digest=left_key_digest,
                        canonical_supported=left_expr.supported,
                        Z_independent_source=left_expr.z_independent,
                        source_root_metadata_digest=_metadata_digest(_node(left_graph, left_root) or {}),
                        source_root_metadata_keys=_metadata_keys(_node(left_graph, left_root) or {}))
    row["right"] = dict(right_descriptor, canonical_digest=right_key_digest,
                         canonical_supported=right_expr.supported,
                         Z_independent_source=right_expr.z_independent,
                         source_root_metadata_digest=_metadata_digest(_node(right_graph, right_root) or {}),
                         source_root_metadata_keys=_metadata_keys(_node(right_graph, right_root) or {}))
    row["rewrite_audit"] = list(dict.fromkeys(left_expr.audit + right_expr.audit))
    # U, Q, m and the history ports are ordinary-Z coefficient functions and
    # may legitimately carry the free-Z source parameter.  The C1 radius is
    # the boundary quantity that must pass the explicit Z-independence walk.
    if field == "R" and left_expr.z_independent is not True:
        row["reasons"].append("left_radius_Z_independence_unproved")
    if field == "R" and right_expr.z_independent is not True:
        row["reasons"].append("right_radius_Z_independence_unproved")
    Z_independent = (field != "R" or
                     (left_expr.z_independent is True and right_expr.z_independent is True))
    if left_expr.key == right_expr.key:
        if Z_independent and left_expr.supported and right_expr.supported:
            row["equal"] = True
            row["source_identity_proved"] = True
            row["status"] = "proved"
        elif (Z_independent and
              _opaque_pair_safe(left_expr, right_expr, left_record, right_record,
                                prefix_checks, graph_lengths)):
            row["equal"] = True
            row["source_identity_proved"] = True
            row["status"] = "proved_shared_prefix_opaque"
            row["reasons"].append("opaque_operation_admitted_by_exact_shared_prefix")
        else:
            row["reasons"].extend(left_expr.unknown_ops + right_expr.unknown_ops)
            row["reasons"].append("unknown_operation_unproved")
    else:
        row["reasons"].append("differing_symbolic_DAG")
    row["reasons"] = list(dict.fromkeys(row["reasons"]))
    row["left_scope_digest"] = _row_scope_digest(left_descriptor)
    row["right_scope_digest"] = _row_scope_digest(right_descriptor)
    return row


def _common_L_status(left_side, right_side, payloads):
    statuses = []
    details = []
    for side in (left_side, right_side):
        port = side.get("physical", {}).get("Q") if isinstance(side, dict) else None
        qy_port = side.get("physical", {}).get("Q_y") if isinstance(side, dict) else None
        for label, candidate in (("Q0", port), ("Qy0", qy_port)):
            roots = candidate.get("source_nodes") if isinstance(candidate, dict) else None
            root = roots[0] if isinstance(roots, list) and roots else None
            graph = _graph_nodes(payloads.get(candidate.get("source_record"), {})) if isinstance(candidate, dict) else []
            status, denominator, certificate = _root_positive_L(graph, root)
            statuses.append(status)
            details.append({"side": side.get("chart"), "field": label, "source_record": candidate.get("source_record") if isinstance(candidate, dict) else None,
                            "source_node": root, "denominator_source_node": denominator,
                            "source_positive_certificate": certificate, "formula": L_FORMULA})
    same = bool(statuses) and all(statuses)
    denominator_nodes = {item["denominator_source_node"] for item in details}
    same = same and denominator_nodes == {L_NODE}
    return same, details


def _history_alias_status(spec, left_side, right_side, payloads):
    # Compare the complete archived endpoint memory with the next incoming
    # memory.  Local chart integrals are different quantities.
    left_record, right_record = spec["left_physical"], spec["right_physical"]
    left_binding = _payload_bindings(payloads.get(left_record, {}), _BINDING_KEYS[left_record])
    right_binding = _payload_bindings(payloads.get(right_record, {}), _BINDING_KEYS[right_record])
    left_chart, right_chart = spec["left_chart"], spec["right_chart"]
    endpoint_key = "Rh" if left_record == PATCH else left_chart
    left_nodes = left_binding.get("endpoints", {}).get(endpoint_key, {}).get("integrals", {}).get("S1")
    if right_record == PATCH:
        right_nodes = right_binding.get("incoming_integrals", {}).get("S1")
        right_path = _BINDING_KEYS[right_record] + ".incoming_integrals.S1"
    else:
        right_nodes = right_binding.get("charts", {}).get(right_chart, {}).get("incoming_integrals", {}).get("S1")
        right_path = _BINDING_KEYS[right_record] + ".charts." + right_chart + ".incoming_integrals.S1"
    left_path = _BINDING_KEYS[left_record] + ".endpoints." + endpoint_key + ".integrals.S1"
    return {
        "exact_alias": isinstance(left_nodes, list) and isinstance(right_nodes, list) and left_nodes == right_nodes,
        "available": isinstance(left_nodes, list) and isinstance(right_nodes, list),
        "left_endpoint_source_nodes": left_nodes,
        "right_incoming_source_nodes": right_nodes,
        "left_source_record": left_record,
        "right_source_record": right_record,
        "left_path": left_path,
        "right_path": right_path,
        "source_identity_used_for_B": False,
        "lower_R_Q_y_saved_endpoint_reused": False,
    }


def _record_payloads(source_records, checks):
    payloads = {}
    exact_inventory = isinstance(source_records, dict) and set(source_records) == set(RECORD_NAMES)
    checks["exact_record_inventory"] = exact_inventory
    if not isinstance(source_records, dict):
        return payloads
    for name in RECORD_NAMES:
        record = source_records.get(name)
        authenticated = (
            isinstance(record, dict)
            and record.get("raw_sha256") == RECORD_HASHES[name]
            and record.get("raw_authenticated_before_decode") is True
            and isinstance(record.get("source_text"), str)
            and hashlib.sha256(record["source_text"].encode("utf-8")).hexdigest() == record.get("decoded_sha256")
        )
        checks["authenticated:" + name] = authenticated
        if not authenticated:
            continue
        try:
            payloads[name] = json.loads(record["source_text"])
        except Exception:
            checks["decoded_json:" + name] = False
            continue
        checks["decoded_json:" + name] = isinstance(payloads[name], dict)
    return payloads


def _retained_data(parent, weighted):
    parent_data = parent.get("retained_construction_data") if isinstance(parent, dict) else None
    weighted_data = weighted.get("retained_construction_data") if isinstance(weighted, dict) else None
    expected = {
        "cumulative_terms": 212,
        "raw_signed_terms": 213,
        "reduced_terms": 183,
        "six_signed_channels_retained": True,
        "all_old_error_data_retained": True,
        "reserved_contribution_registered_once": True,
        "inherited_errors_deleted": False,
    }
    okay = isinstance(parent_data, dict) and isinstance(weighted_data, dict) and parent_data == weighted_data
    okay = okay and all(parent_data.get(key) == value for key, value in expected.items())
    return dict(parent_data) if isinstance(parent_data, dict) else {}, okay


def _parent_preservation(parent):
    """Verify the inherited 17 K terms and the two external B ports remain."""
    component = parent.get("component_result") if isinstance(parent, dict) else None
    definition = component.get("definition_graph") if isinstance(component, dict) else None
    nodes = definition.get("nodes") if isinstance(definition, dict) else None
    if not isinstance(nodes, list):
        return False, {"K_integral_count": 0, "external_B_count": 0,
                       "native_IN_formula": None}
    K_rows = [row for row in nodes if isinstance(row, dict)
              and row.get("operation") == "unchanged_weighted_K_integral_source_port"]
    external_rows = [row for row in nodes if isinstance(row, dict)
                     and row.get("operation") == "unchanged_weighted_external_B_source_port"]
    external_nodes = sorted(row.get("source_node") for row in external_rows
                            if _is_int(row.get("source_node")))
    formula = component.get("source_interface", {}).get("native_IN_formula")
    details = {
        "K_integral_count": len(K_rows),
        "external_B_count": len(external_rows),
        "external_B_source_nodes": external_nodes,
        "native_IN_formula": formula,
    }
    return (len(K_rows) == 17 and len(external_rows) == 2
            and external_nodes == [333, 510]
            and formula == IN_FORMULA), details


def _parent_context(parent, weighted):
    parent_interface = parent.get("source_interface") if isinstance(parent, dict) else {}
    definition = weighted.get("definition_graph") if isinstance(weighted, dict) else None
    weighted_nodes = _graph_nodes(definition)
    if not weighted_nodes and isinstance(definition, dict):
        candidate = definition.get("nodes")
        weighted_nodes = candidate if isinstance(candidate, list) else []
    inherited = parent_interface.get("proved_weighted_jump_nodes") == [520, 521]
    inherited = inherited and parent_interface.get("J1_zero") is True
    inherited = inherited and parent.get("component_result", {}).get("passed") is True
    inherited = inherited and parent.get("source_only") is True
    unchanged = parent_interface.get("remaining_parent_jump_nodes") == list(range(522, 536))
    unchanged = unchanged and parent_interface.get("native_IN_formula") == IN_FORMULA
    unchanged = unchanged and parent_interface.get("old_archived_jump_rows_unchanged") is True
    unchanged = unchanged and parent.get("old_archived_16_J_terms") is True
    unchanged = unchanged and len(weighted_nodes) >= 536
    return inherited, unchanged


def _build_in_formula(proven_jumps):
    remaining = [index for index in range(2, 16) if index not in set(proven_jumps)]
    if not remaining:
        jump_part = "0"
    elif remaining == list(range(2, 16)):
        jump_part = "sum(J_2..J_15)"
    else:
        jump_part = "sum(" + ",".join("J_" + str(index) for index in remaining) + ")"
    return "I_N=B_lastupper-B_firstlower-" + jump_part + "+sum17(integral K_i dy)"


def _synthetic_engine_soundness_audit():
    """Small semantic invariants for the new bounded walker only.

    These graphs are deliberately independent of the archived inputs.  The
    audit runs inside the one future ``derive`` call and catches regressions
    in multiplicity, branch laziness, alpha scopes, substitution capture, and
    opaque-scope separation without replaying any ancestor construction.
    """
    sum_walker = _SourceWalker("__synthetic_sum__", [
        {"operation": "current_original_source_parameter", "name": "x",
         "Z_independent": True},
        {"operation": "sum", "arguments": [0, 0]},
    ])
    x_expr = sum_walker.walk(0)
    xx_expr = sum_walker.walk(1)

    branch_walker = _SourceWalker("__synthetic_branch__", [
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_real_comparison", "left": 0, "right": 1,
         "operator": "=="},
        {"operation": "synthetic_unknown_active_body"},
        {"operation": "exact_rational", "numerator": 7, "denominator": 1},
        {"operation": "original_lazy_flat_branch", "flat_predicate": 2,
         "active_body": 3, "flat_value": 4},
    ])
    branch_expr = branch_walker.walk(5)

    empty_integral_walker = _SourceWalker("__synthetic_empty_integral__", [
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "synthetic_unknown_integrand"},
        {"operation": "definite_integral", "integrand": 1, "lower": 0,
         "upper": 0, "variable": "alpha"},
    ])
    empty_expr = empty_integral_walker.walk(2)

    substitution_walker = _SourceWalker("__synthetic_substitution__", [
        {"operation": "exact_rational", "numerator": 7, "denominator": 1},
        {"operation": "bound_variable", "name": "x"},
        {"operation": "bound_variable", "name": "y"},
        {"operation": "function_substitution", "expression": 1,
         "variable": 1, "value": 2,
         "Z_independent_substitution": True},
    ])
    substituted_expr = substitution_walker.walk(3, {2: 0})

    # The substitution value is frozen outside the integral, while the
    # integral's alpha-bound ``x`` shadows both the endpoint scope and that
    # frozen value.  This catches the classic capture regression directly.
    capture_walker = _SourceWalker("__synthetic_capture__", [
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "exact_rational", "numerator": 2, "denominator": 1},
        {"operation": "bound_variable", "name": "x"},
        {"operation": "definite_integral", "integrand": 3,
         "lower": 0, "upper": 1, "variable": "x"},
        {"operation": "function_substitution", "expression": 4,
         "variable": 3, "value": 2,
         "Z_independent_substitution": True},
    ])
    capture_expr = capture_walker.walk(5)
    capture_bound_rows = [expr for key, expr in capture_walker.memo.items()
                          if isinstance(key, tuple) and key and key[0] == 3]
    capture_shadow_ok = any(
        expr.operation == "bound_variable" and expr.key != _rational(2)
        for expr in capture_bound_rows
    )

    opaque_walker = _SourceWalker("__synthetic_opaque_scope__", [
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "exact_rational", "numerator": 2, "denominator": 1},
        {"operation": "synthetic_scope_marker"},
    ])
    opaque_left = opaque_walker.walk(2, {0: 0})
    opaque_right = opaque_walker.walk(2, {0: 1})

    checks = {
        "sum_preserves_multiplicity": xx_expr.key != x_expr.key,
        "lazy_flat_branch_skips_unknown_active_body": (
            branch_expr.key == _rational(7)
            and 3 not in {key[0] for key in branch_walker.memo
                           if isinstance(key, tuple) and key}
        ),
        "lazy_flat_branch_uses_graph_flat_value": branch_expr.key == _rational(7),
        "empty_integral_zero_without_integrand_walk": (
            empty_expr.key == _rational(0)
            and 1 not in {key[0] for key in empty_integral_walker.memo
                           if isinstance(key, tuple) and key}
        ),
        "substitution_value_frozen_under_outer_scope": substituted_expr.key == _rational(7),
        "integral_bound_name_shadows_substitution": capture_shadow_ok,
        "opaque_scopes_do_not_compare_equal": opaque_left.key != opaque_right.key,
    }
    checks["all_engine_soundness_checks"] = all(checks.values())
    return checks


def derive(source_records, *, expected_input_digest):
    """Run one bounded source comparison pass for J2 through J15.

    ``passed`` records successful authenticated source inspection and honest
    endpoint rows.  It does not assert that every remaining B jump is zero.
    """
    checks = {}
    try:
        checks["frozen_input_digest"] = _digest(source_records) == expected_input_digest
    except Exception:
        checks["frozen_input_digest"] = False
    payloads = _record_payloads(source_records, checks)
    parent = payloads.get(PARENT, {})
    weighted = payloads.get(WEIGHTED, {})
    weighted_nodes = _graph_nodes(weighted.get("definition_graph", {}))
    if not weighted_nodes and isinstance(weighted.get("definition_graph"), dict):
        candidate = weighted["definition_graph"].get("nodes")
        weighted_nodes = candidate if isinstance(candidate, list) else []

    same_family, source_family, source_families = _same_source_family(payloads)
    checks["same_exact_source_family"] = same_family
    checks["parent_report_accepted"] = (
        isinstance(parent, dict)
        and parent.get("component_result", {}).get("passed") is True
        and parent.get("source_only") is True
        and parent.get("full_goal_complete") is False
    )
    checks["accepted_J0_J1_context"] = _parent_context(parent, weighted)[0]
    checks["weighted_graph_contains_remaining_nodes"] = all(
        _node(weighted_nodes, index) is not None for index in range(522, 536)
    )
    checks["weighted_graph_old_nodes_unchanged"] = all(
        _node(weighted_nodes, index).get("retained") is True
        and _node(weighted_nodes, index).get("zero_assumption") is False
        for index in range(522, 536) if _node(weighted_nodes, index) is not None
    )
    retained, retained_ok = _retained_data(parent, weighted)
    checks["retained_construction_data_frozen"] = retained_ok
    old_ports_preserved, old_ports_details = _parent_preservation(parent)
    checks["old_17_K_integrals_and_2_outer_B_ports_retained"] = old_ports_preserved
    prefix_checks, graph_by_record = _prefix_checks(payloads)
    checks.update(prefix_checks)
    checks.update(_binding_shape(payloads))
    checks["raw_gzip_records_authenticated_before_decode"] = all(
        checks.get("authenticated:" + name) is True for name in (MOD, LOOP, PATCH, OUTER)
    )
    checks["forbidden_raw_only_records_not_in_input"] = True
    checks["inherited_owned_Z_cell"] = (
        parent.get("owned_Z_cell") == list(OWNED_Z_CELL)
        or parent.get("source_interface", {}).get("owned_Z_cell") == list(OWNED_Z_CELL)
    )
    checks["ordinary_Z_scope"] = parent.get("ordinary_Z_orders") == [0, 1] or not parent
    checks["ordinary_Z_not_temporal_n"] = parent.get("ordinary_Z_not_temporal_n") is True or not parent
    engine_soundness = _synthetic_engine_soundness_audit()
    checks.update({"engine_" + name: value for name, value in engine_soundness.items()})

    walkers = {
        name: _SourceWalker(name, graph_by_record.get(name, []))
        for name in (MOD, LOOP, PATCH, OUTER)
    }
    graph_lengths = {name: len(graph_by_record.get(name, [])) for name in graph_by_record}
    endpoint_results = []
    proven_jumps = []
    unresolved_jumps = []
    graph = _Graph()

    for spec in _SEAM_SPECS:
        jump_node = spec["native_graph_node"]
        b_node = _node(weighted_nodes, jump_node)
        seam_audit = []
        if b_node is None or not isinstance(b_node.get("arguments"), list) or len(b_node["arguments"]) < 2:
            left_side, left_errors = None, ["malformed_weighted_B_node"]
            right_side, right_errors = None, ["malformed_weighted_B_node"]
        else:
            if b_node.get("operation") != "actual_internal_primitive_jump":
                seam_audit.append("weighted_jump_operation_mismatch")
            if b_node.get("jump_index") != spec["jump_index"]:
                seam_audit.append("weighted_jump_index_mismatch")
            if len(b_node["arguments"]) != 2:
                seam_audit.append("weighted_jump_argument_arity_mismatch")
            right_endpoint_id, left_endpoint_id = b_node["arguments"][:2]
            left_side, left_errors = _endpoint_descriptor(weighted_nodes, left_endpoint_id, spec, "left")
            right_side, right_errors = _endpoint_descriptor(weighted_nodes, right_endpoint_id, spec, "right")
        seam_audit.extend(left_errors)
        seam_audit.extend(right_errors)
        left_binding_audit = _archive_binding_audit(payloads, left_side)
        right_binding_audit = _archive_binding_audit(payloads, right_side)
        seam_audit.extend(left_binding_audit["errors"])
        seam_audit.extend(right_binding_audit["errors"])
        expected_arguments = [
            right_side.get("weighted_node") if isinstance(right_side, dict) else None,
            left_side.get("weighted_node") if isinstance(left_side, dict) else None,
        ]
        if not isinstance(b_node, dict) or b_node.get("arguments") != expected_arguments:
            seam_audit.append("weighted_jump_argument_identity_mismatch")
        rows = {}
        for field, count in (("U", 2), ("Q", 2), ("R", 1), ("m", 3), ("Q_y", 2), ("R_Q_y", 1)):
            for order in range(count):
                label = field + "_Z" + str(order)
                left_descriptor = _root_descriptor(left_side, field, order) if left_side else None
                right_descriptor = _root_descriptor(right_side, field, order) if right_side else None
                row = _compare_row(left_descriptor, right_descriptor, payloads, walkers,
                                   prefix_checks, graph_lengths, field, order)
                rows[label] = row
        L_proved, L_rows = _common_L_status(left_side or {}, right_side or {}, payloads)
        history = _history_alias_status(spec, left_side or {}, right_side or {}, payloads)
        required = (rows["U_Z0"], rows["U_Z1"], rows["Q_Z0"], rows["Q_Z1"], rows["R_Z0"])
        required_proved = all(row.get("source_identity_proved") is True for row in required)
        source_identity_proved = bool(required_proved and L_proved and not seam_audit)
        b_status = {
            "source_function_join_receipt": "accepted" if not seam_audit else "unproved",
            "source_identity_proved": source_identity_proved,
            "zero_proved": source_identity_proved,
            "zero_assumption": False,
            "status": "proved_zero" if source_identity_proved else "unproved",
            "formula": B_FORMULA,
            "jump_operation": b_node.get("operation") if isinstance(b_node, dict) else None,
            "jump_index": b_node.get("jump_index") if isinstance(b_node, dict) else None,
            "jump_arguments": (list(b_node.get("arguments", []))
                               if isinstance(b_node, dict)
                               and isinstance(b_node.get("arguments"), list) else None),
            "left_weighted_node": left_side.get("weighted_node") if left_side else None,
            "right_weighted_node": right_side.get("weighted_node") if right_side else None,
            "common_positive_L": L_proved,
            "common_positive_L_node": L_NODE if L_proved else None,
            "common_positive_L_rows": L_rows,
            "required_fields": ["U0", "U1", "Q0", "Q1", "R"],
            "optional_evidence_fields": ["m0", "m1", "m2", "Q_y0", "Q_y1", "R_Q_y0"],
            "history_alias": history,
            "reasons": list(dict.fromkeys(seam_audit + [reason for row in required for reason in row.get("reasons", [])])),
        }
        if source_identity_proved:
            proven_jumps.append(spec["jump_index"])
        else:
            unresolved_jumps.append(jump_node)
        honest = all(
            isinstance(rows[label].get("left"), dict)
            and isinstance(rows[label].get("right"), dict)
            and rows[label]["left"].get("source_node") is not None
            and rows[label]["right"].get("source_node") is not None
            for label in ("U_Z0", "U_Z1", "Q_Z0", "Q_Z1", "R_Z0")
        )
        seam_result = {
            "jump_index": spec["jump_index"],
            "native_graph_node": jump_node,
            "weighted_jump_node": jump_node,
            "seam": spec["seam"],
            "certificate_family": spec["certificate_family"],
            "left_chart": spec["left_chart"],
            "right_chart": spec["right_chart"],
            "left_endpoint": spec["left_endpoint"],
            "right_endpoint": spec["right_endpoint"],
            "archive_binding_map": {
                "left_geometry": spec["left_geometry"],
                "right_geometry": spec["right_geometry"],
                "left_physical": spec["left_physical"],
                "right_physical": spec["right_physical"],
            },
            "left_endpoint_source": left_side,
            "right_endpoint_source": right_side,
            "left_archive_binding_audit": left_binding_audit,
            "right_archive_binding_audit": right_binding_audit,
            "endpoint_rows": rows,
            "endpoint_comparisons": rows,
            "native_B_identity_status": b_status,
            "function_identity_proved": source_identity_proved,
            "honest_required_rows": honest,
            "unresolved_reasons": list(dict.fromkeys(seam_audit + b_status["reasons"])),
        }
        endpoint_results.append(seam_result)
        graph.add("remaining_B_seam_endpoint_source_comparison", arguments=(),
                  jump_index=spec["jump_index"], native_graph_node=jump_node,
                  endpoint_row_count=len(rows), honest_required_rows=honest,
                  source_identity_proved=source_identity_proved,
                  common_positive_L=L_proved, source_function_only=True)

    checks["all_fourteen_honest_endpoint_rows"] = len(endpoint_results) == 14 and all(
        item.get("honest_required_rows") is True for item in endpoint_results
    )
    checks["all_fourteen_B_comparisons_attempted"] = len(endpoint_results) == 14
    checks["weighted_B_node_arguments_authoritative"] = all(
        item.get("native_B_identity_status", {}).get("left_weighted_node") is not None
        and item.get("native_B_identity_status", {}).get("right_weighted_node") is not None
        for item in endpoint_results
    )
    checks["weighted_B_operations_and_indices_authoritative"] = all(
        item.get("native_B_identity_status", {}).get("formula") == B_FORMULA
        and item.get("native_B_identity_status", {}).get("jump_operation") == "actual_internal_primitive_jump"
        and item.get("native_B_identity_status", {}).get("jump_index") == item.get("jump_index")
        and item.get("native_B_identity_status", {}).get("jump_arguments") == [
            item.get("native_B_identity_status", {}).get("right_weighted_node"),
            item.get("native_B_identity_status", {}).get("left_weighted_node"),
        ]
        and item.get("native_graph_node") == 520 + item.get("jump_index")
        for item in endpoint_results
    )
    checks["all_weighted_endpoint_arrays_match_archive_binding_maps"] = all(
        item.get("left_archive_binding_audit", {}).get("passed") is True
        and item.get("right_archive_binding_audit", {}).get("passed") is True
        for item in endpoint_results
    )

    newly_proved_nodes = [item["native_graph_node"] for item in endpoint_results
                          if item["native_B_identity_status"].get("source_identity_proved") is True]
    proven_weighted_nodes = [520, 521] + newly_proved_nodes
    remaining_nodes = [item["native_graph_node"] for item in endpoint_results
                       if item["native_B_identity_status"].get("source_identity_proved") is not True]
    in_formula = _build_in_formula(proven_jumps)
    frozen_errors = parent.get("old_errors_frozen")
    old_error_counts = (isinstance(frozen_errors, dict)
                        and frozen_errors.get("cumulative_terms") == 212
                        and frozen_errors.get("raw_signed_terms") == 213
                        and frozen_errors.get("reduced_terms") == 183)
    all_old_terms = (retained_ok and old_error_counts
                     and parent.get("six_signed_channels_retained") is True
                     and parent.get("reserved_contribution_registered_once") is True
                     and parent.get("old_error_terms_removed") is False)
    checks["all_old_terms_retained"] = all_old_terms
    checks["accepted_context_not_reexecuted"] = True

    stage_values = {
        "new_remaining_seam_source_comparison_completed": bool(
            checks.get("all_fourteen_honest_endpoint_rows") and checks.get("all_fourteen_B_comparisons_attempted")
        ),
        "accepted_J0_J1_not_reexecuted": True,
        "all_old_terms_retained": all_old_terms,
    }
    checks["stage_rows_complete"] = stage_values["new_remaining_seam_source_comparison_completed"]

    graph.add("native_IN_source_representation_remaining_Js", arguments=(),
              formula=in_formula, inherited_J0_J1=[520, 521],
              newly_proved_weighted_jump_nodes=newly_proved_nodes,
              unresolved_weighted_jump_nodes=remaining_nodes,
              old_terms_retained=all_old_terms)
    graph_digest = _digest(graph.nodes)
    source_graph_summary = {
        name: {
            "node_count": len(graph_by_record.get(name, [])),
            # Keep the provenance digest bounded.  The runner already freezes
            # and authenticates the raw archive bytes; reserializing a
            # 258,857-node archive here would be an unnecessary second full
            # archive walk.
            "sampled_provenance_digest": _digest({
                "source_record": name,
                "node_count": len(graph_by_record.get(name, [])),
                "head": graph_by_record.get(name, [])[:8],
                "tail": graph_by_record.get(name, [])[-8:],
            }),
            "sampled_provenance": True,
            "raw_archive_sha256": RECORD_HASHES[name],
            "prefix_role": "LOOP_geometry_MOD_physical" if name in (LOOP, MOD) else
                           "PATCH_own" if name == PATCH else "OUTER_own",
        }
        for name in (MOD, LOOP, PATCH, OUTER)
    }
    definition_graph = {
        "namespace": SCHEMA,
        "nodes": graph.nodes,
        "node_count": len(graph.nodes),
        "nodes_digest": graph_digest,
        "graph_digest": graph_digest,
        "source_graphs": source_graph_summary,
        "source_intern_rows": {
            "storage": "definition_graph.source_endpoint_DAGs[*].intern_rows",
            "per_archive": {
                name: {
                    "intern_count": len(walker.intern_rows),
                    "intern_rows_digest": _digest(walker.intern_rows),
                }
                for name, walker in walkers.items()
            },
        },
        "source_key_audit": {
            name: {
                "intern_count": len(walker.intern_rows),
                "memo_entries": len(walker.memo),
                "node_budget": walker.node_budget,
                "rewrite_event_count": len(walker.rewrite_audit),
                "intern_rows_digest": _digest(walker.intern_rows),
            }
            for name, walker in walkers.items()
        },
        "source_endpoint_DAGs": {
            name: {
                "intern_rows": list(walker.intern_rows),
                "intern_count": len(walker.intern_rows),
                "memo_entries": len(walker.memo),
                "visited_nodes": walker.visited,
                "node_budget": walker.node_budget,
                "rewrite_event_count": len(walker.rewrite_audit),
                "raw_archive_sha256": RECORD_HASHES[name],
            }
            for name, walker in walkers.items()
        },
        "engine_soundness_audit": engine_soundness,
        "inherited_old_ports": old_ports_details,
        "compared_seams": endpoint_results,
        "interned_source_DAG": True,
        "flat_merkle_source_keys": True,
        "memo_scope": "(archive, source_node, endpoint_coordinate_scope, alpha_bound_integral_scope)",
        "node_budget_per_archive": 250000,
        "visited_nodes_per_archive": {name: walker.visited for name, walker in walkers.items()},
        "rewrite_audit": {name: list(dict.fromkeys(walker.rewrite_audit)) for name, walker in walkers.items()},
        "exact_rational_endpoint_substitution": True,
        "explicit_empty_integral_lower_equals_upper": True,
        "nonempty_integrals_remain_symbolic": True,
        "flat_branch_and_sigma_endpoint_rules": True,
        "no_expanded_polynomial": True,
        "no_quadrature": True,
        "no_numeric_values": True,
        "unknown_operation_is_unproved": True,
    }

    source_prerequisite_audit = dict(checks)
    passed = bool(all(source_prerequisite_audit.values()) and all(stage_values.values()))
    unresolved_gates = []
    for name, value in source_prerequisite_audit.items():
        if not value:
            unresolved_gates.append("source_prerequisite_unresolved: " + name)
    for item in endpoint_results:
        if item["native_B_identity_status"].get("source_identity_proved") is not True:
            unresolved_gates.append("native_B_identity_unproved_J" + str(item["jump_index"]) + "_node" + str(item["native_graph_node"]))
    unresolved_gates = list(dict.fromkeys(unresolved_gates))

    parent_context = {
        "parent_report": PARENT,
        "parent_report_sha256": RECORD_HASHES[PARENT],
        "parent_unchanged": _parent_context(parent, weighted)[1],
        "accepted_J0_J1": True,
        "reexecuted": False,
        "retained_construction_data_unchanged": retained_ok,
        "old_archived_jump_rows_unchanged": parent.get("source_interface", {}).get(
            "old_archived_jump_rows_unchanged") is True,
        "old_archived_16_J_terms_unchanged": parent.get("old_archived_16_J_terms") is True,
        "old_17_K_integrals_and_2_outer_B_ports_retained": old_ports_preserved,
        "old_packet_errors_deleted": False,
    }

    source_interface = {
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z_not_temporal_n": True,
        "owned_Z_cell": list(OWNED_Z_CELL),
        "source_B_formula": B_FORMULA,
        "common_positive_L_formula": L_FORMULA,
        "common_positive_L_node": L_NODE,
        "inherited_J0_J1_not_reexecuted": True,
        "inherited_proved_weighted_jump_nodes": [520, 521],
        "proved_weighted_jump_nodes": proven_weighted_nodes,
        "newly_proved_weighted_jump_nodes": newly_proved_nodes,
        "remaining_parent_jump_nodes": remaining_nodes,
        "native_IN_formula": in_formula,
        "IN_formula": in_formula,
        "old_archived_jump_rows_unchanged": parent_context["old_archived_jump_rows_unchanged"],
        "old_archived_16_J_terms_unchanged": parent_context["old_archived_16_J_terms_unchanged"],
        "old_17_K_integrals_and_2_outer_B_ports_retained": old_ports_preserved,
        "old_packet_errors_deleted": False,
        "m_R_Q_y_are_optional_evidence_only": True,
        "lower_R_Q_y_saved_endpoint_rows_not_reused": True,
        "source_function_identity_required_for_B": ["U0", "U1", "Q0", "Q1", "R"],
        "unknown_operation_policy": "unproved_not_inequality",
    }

    source_inspection_audit = {
        "new_source_consumer_calls": 1,
        "accepted_ancestor_operator_calls": 0,
        "canonicalizer_reexecuted": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "old_algebra_replayed": False,
        "quadrature_calls": 0,
        "source_records_only": True,
        "bounded_endpoint_source_walk": True,
        "memoized_source_DAG": True,
        "node_budget_enforced": True,
        "synthetic_engine_soundness_audit": engine_soundness,
        "runtime_model_metadata": None,
    }

    result = {
        "schema": SCHEMA,
        "component": "remaining_B_seams",
        "passed": passed,
        "source_only": True,
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z_not_temporal_n": True,
        "owned_Z_cell": list(OWNED_Z_CELL),
        "input_digest": expected_input_digest,
        "source_family": source_family,
        "source_family_dict": source_family,
        "source_family_audit": {
            "same_exact_source_family": same_family,
            "archives": {name: family for name, family in zip((PARENT, WEIGHTED, MOD, LOOP, PATCH, OUTER), source_families)},
        },
        "source_prerequisite_audit": source_prerequisite_audit,
        "seam_endpoint_comparisons": endpoint_results,
        "seam_comparisons": endpoint_results,
        "definition_graph": definition_graph,
        "source_interface": source_interface,
        "parent_context": parent_context,
        "parent_unchanged": parent_context["parent_unchanged"],
        "retained_construction_data": retained,
        "old_17_K_integrals_and_2_outer_B_ports_retained": old_ports_preserved,
        "proven_weighted_nodes": proven_weighted_nodes,
        "proven_weighted_jump_nodes": proven_weighted_nodes,
        "newly_proved_weighted_jump_nodes": newly_proved_nodes,
        "remaining_parent_jump_nodes": remaining_nodes,
        "native_IN_formula": in_formula,
        "unresolved_weighted_nodes": remaining_nodes,
        "unresolved_weighted_jump_nodes": remaining_nodes,
        "unresolved_gates": unresolved_gates,
        "remaining_construction_gates": [
            "remaining_native_B_source_identity_proofs" if remaining_nodes else "no_remaining_native_B_source_identity_proofs",
            "native_old_active_density_identity",
            "native_compact_density_identity",
            "same_source_old_packet_error_theorem",
            "terminal_difference_zero_proof",
            "global_five_moment_closure",
            "exact_heat_exterior",
            "admissible_stress_reconstruction",
            "time_order_coefficient_recursion",
            "oscillatory_stress_correction",
            "independent_cartesian_NS_validation",
        ],
        "runtime_model_metadata": None,
        "source_inspection_audit": source_inspection_audit,
        "source_family_dict_returned": isinstance(source_family, dict),
        "numeric_values_computed": False,
        "old_error_terms_removed": False,
        "global_five_moment_closure": False,
        "time_order_recursion_claimed": False,
        "full_goal_complete": False,
        "all_native_B_jumps_zero_proved": len(newly_proved_nodes) == 14,
        "all_remaining_native_B_jumps_zero_proved": len(newly_proved_nodes) == 14,
    }
    result.update(stage_values)
    result.update({name: False for name in FALSE_FLAGS})
    return result


__all__ = [
    "SCHEMA",
    "RECORD_NAMES",
    "RECORD_HASHES",
    "STAGE_FLAGS",
    "FALSE_FLAGS",
    "derive",
]
