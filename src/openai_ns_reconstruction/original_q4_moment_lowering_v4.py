"""V4 source-preserving, hygienic lowering for the archived original Q4 moment.

This module deliberately implements an IR compiler, not an evaluator.  Graph
rows and graph handles are retained in the result so that a later consumer can
choose the runtime semantics.  In particular, opaque source rows are never
turned into constants and exact rationals are never algebraically folded.

The public entry point is :func:`build`.  It accepts a small graph adapter
shape (a mapping containing ``nodes`` or a direct handle-to-row mapping) and
returns JSON-compatible dictionaries only.
"""

# The compiler is intentionally self-contained.  It uses no project imports:
# this keeps the source archive independent from any evaluator implementation.

_MAX_NODES = 512
_MAX_DEPTH = 128
_MAX_ENV = 128

_OUTER_VARIABLE = "original_n1_true_core_rho"
_INTEGRAL_VARIABLE = "original_core_inlet_moment_t"
_EXPECTED_FAMILY_NAMES = {
    "q4",
    "original_q4",
    "q4_moment",
    "original_q4_moment",
    "original-q4",
    "original-q4-moment",
}

_STANDARD_OPS = {
    "sum",
    "product",
    "negative",
    "positive_quotient",
    "analytic_unary",
    "exact_rational",
    "bound_variable",
    "function_substitution",
    "definite_integral",
}

_OPAQUE_OPS = {
    "current_original_source_parameter",
    "current_original_unique_analytic_core_component",
}

_ANALYTIC_UNARY_NAMES = {
    "abs",
    "cos",
    "cosh",
    "erf",
    "exp",
    "expm1",
    "log",
    "log1p",
    "sin",
    "sinh",
    "sqrt",
    "tan",
    "tanh",
}

_MISSING = object()


def _clone(value, depth=0):
    """Clone JSON-like input without importing a serializer or evaluator."""

    if depth > _MAX_DEPTH:
        raise ValueError("clone depth limit")
    if isinstance(value, dict):
        return {key: _clone(item, depth + 1) for key, item in value.items()}
    if isinstance(value, list):
        return [_clone(item, depth + 1) for item in value]
    if isinstance(value, tuple):
        return [_clone(item, depth + 1) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    # A graph handle is allowed to be a scalar object only when it is already
    # JSON-like.  Rejecting other objects makes the output JSON-ready.
    raise ValueError("non-json graph value")


def _first(row, names, default=_MISSING):
    for name in names:
        if isinstance(row, dict) and name in row:
            return row[name]
    return default


def _row_op(row):
    if not isinstance(row, dict):
        return None
    op = _first(row, ("operation", "op", "operator", "kind", "type"), None)
    return op if isinstance(op, str) else None


def _ref_key(handle):
    """Return the JSON object key used by the source-atom report."""

    if isinstance(handle, bool):
        return "true" if handle else "false"
    if isinstance(handle, str):
        return handle
    if isinstance(handle, (int, float)) or handle is None:
        return str(handle)
    return None


def _family_name(family):
    if isinstance(family, str):
        return family.strip().lower().replace(" ", "_")
    if isinstance(family, dict):
        label = _first(family, ("name", "family", "source_family", "label", "id"), _MISSING)
        if isinstance(label, str):
            return label.strip().lower().replace(" ", "_")
    return None


def _family_sha256_values(family):
    if not isinstance(family, dict):
        return []
    values = []
    for key, value in family.items():
        key_text = str(key).lower()
        if "sha256" in key_text and isinstance(value, str):
            values.append(value)
        elif isinstance(value, str) and len(value) == 64:
            # Archive family records sometimes use concise field names while
            # retaining the SHA-256 values themselves.
            values.append(value)
    return values


def _family_is_valid(family):
    if isinstance(family, dict):
        for key in ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256'):
            value = family.get(key)
            if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdefABCDEF' for c in value):
                return False
        label = _family_name(family)
        if label is not None and ("p1" in label or "m15" in label):
            return False
        return True
    return False


def _family_equal(left, right):
    return isinstance(left, dict) and isinstance(right, dict) and left == right


class _Failure(Exception):
    """Internal fail-closed signal carrying a stable, JSON-safe reason."""

    def __init__(self, reason, detail=None):
        self.reason = reason
        self.detail = detail
        Exception.__init__(self, reason)


class _Compiler:
    def __init__(self, graph, inlet_handle, family):
        self.graph = graph
        self.inlet_handle = inlet_handle
        self.family = family
        self.family_name = _family_name(family)
        self.nodes = self._extract_nodes(graph)
        self.metadata = self._extract_metadata(graph)
        self.active = []
        self.seen = 0
        self.binder_count = 0
        self.source_atoms = []
        self.source_atom_map = {}
        self.reasons = []
        self.root_handle = inlet_handle
        self.outer_handle = None
        self.integral_handles = []
        self.scope_candidate_handles = []

    def _extract_nodes(self, graph):
        if not isinstance(graph, dict):
            raise _Failure("invalid_graph")
        candidate = graph.get("exact_graph_nodes", _MISSING)
        if candidate is _MISSING:
            candidate = graph.get("nodes", _MISSING)
        if candidate is _MISSING:
            candidate = graph.get("graph", _MISSING)
            if isinstance(candidate, dict) and "nodes" in candidate:
                candidate = candidate["nodes"]
        if candidate is _MISSING:
            # A direct handle-to-row mapping is useful for tiny archived graphs.
            candidate = {
                key: value
                for key, value in graph.items()
                if isinstance(value, dict) and _row_op(value) is not None
            }
            if not candidate:
                candidate = graph
        if isinstance(candidate, dict):
            return candidate
        if isinstance(candidate, list):
            if candidate and all(
                isinstance(row, dict)
                and _first(row, ("handle", "id", "node_id", "graph_id", "source_ref", "ref", "key"), _MISSING)
                is _MISSING
                for row in candidate
            ):
                # The archived exact_graph_nodes representation is an indexed
                # row list.  A numeric graph handle addresses the list
                # directly; no string-to-node identity is inferred.
                return candidate
            result = {}
            for row in candidate:
                if not isinstance(row, dict):
                    raise _Failure("invalid_node_row")
                handle = _first(
                    row,
                    ("handle", "id", "node_id", "graph_id", "source_ref", "ref", "key"),
                    _MISSING,
                )
                if handle is _MISSING:
                    raise _Failure("node_row_missing_handle")
                key = _ref_key(handle)
                if key is None:
                    raise _Failure("non_json_handle")
                if key in result:
                    raise _Failure("duplicate_node_handle")
                result[handle] = row
            return result
        raise _Failure("invalid_nodes")

    def _extract_metadata(self, graph):
        if not isinstance(graph, dict):
            return {}
        metadata = graph.get("metadata", {})
        if isinstance(metadata, dict):
            result = metadata
        else:
            result = {}
        # Top-level archive metadata remains authoritative when present.
        for name in (
            "source_family",
            "family",
            "bounds",
            "old_errors",
            "errors",
            "numeric",
            "numeric_field",
            "numeric_field_value",
        ):
            if name in graph:
                result = dict(result)
                result[name] = graph[name]
        return result

    def _lookup(self, handle):
        if isinstance(self.nodes, list):
            if isinstance(handle, int) and not isinstance(handle, bool) and 0 <= handle < len(self.nodes):
                return self.nodes[handle]
            raise _Failure("missing_node", {"handle": _clone(handle)})
        if isinstance(self.nodes, dict):
            if handle in self.nodes:
                return self.nodes[handle]
            key = _ref_key(handle)
            if key is not None and key in self.nodes:
                return self.nodes[key]
            # Some JSON graph stores stringified integer keys.  This fallback
            # is intentionally limited to scalar handles.
            if isinstance(handle, (str, int)) and str(handle) in self.nodes:
                return self.nodes[str(handle)]
        raise _Failure("missing_node", {"handle": _clone(handle)})

    def _source_ref(self, handle):
        try:
            return _clone(handle)
        except ValueError:
            raise _Failure("non_json_handle")

    def _record_atom(self, handle, row, op):
        ref = self._source_ref(handle)
        key = _ref_key(ref)
        if key is None:
            raise _Failure("non_json_handle")
        if key not in self.source_atom_map:
            atom = {
                "source_ref": ref,
                "op": op,
                "row": _clone(row),
            }
            self.source_atoms.append(atom)
            self.source_atom_map[key] = atom

    def _node(self, handle):
        self.seen += 1
        if self.seen > _MAX_NODES:
            raise _Failure("node_limit_exceeded")
        if handle in self.active:
            raise _Failure("cycle_detected", {"handle": self._source_ref(handle)})
        row = self._lookup(handle)
        if not isinstance(row, dict):
            raise _Failure("invalid_node_row", {"handle": self._source_ref(handle)})
        try:
            _clone(row)
        except ValueError:
            raise _Failure("non_json_row", {"handle": self._source_ref(handle)})
        return row

    def _push(self, handle):
        if len(self.active) >= _MAX_DEPTH:
            raise _Failure("depth_limit_exceeded")
        if handle in self.active:
            raise _Failure("cycle_detected", {"handle": self._source_ref(handle)})
        self.active.append(handle)

    def _pop(self):
        if self.active:
            self.active.pop()

    def _compile_ref(self, value, env, depth):
        """Compile a graph handle, or a literal inline row when supplied."""

        if depth > _MAX_DEPTH:
            raise _Failure("depth_limit_exceeded")
        if isinstance(value, dict):
            # Inline rows are accepted only as explicit rows, never as arbitrary
            # literal expressions.  This preserves source information.
            if _row_op(value) is None:
                raise _Failure("invalid_inline_reference")
            return self._compile_row(None, value, env, depth + 1, inline=True)
        return self._compile_handle(value, env, depth + 1)

    def _compile_handle(self, handle, env, depth):
        row = self._node(handle)
        self._push(handle)
        try:
            return self._compile_row(handle, row, env, depth, inline=False)
        finally:
            self._pop()

    def _base(self, op, handle, row):
        result = {"op": op}
        if handle is not None:
            result["source_ref"] = self._source_ref(handle)
        result["source_row"] = _clone(row)
        return result

    def _compile_row(self, handle, row, env, depth, inline=False):
        if depth > _MAX_DEPTH:
            raise _Failure("depth_limit_exceeded")
        op = _row_op(row)
        if op is None:
            raise _Failure("invalid_row", {"handle": self._source_ref(handle) if handle is not None else None})
        if op not in _STANDARD_OPS and op not in _OPAQUE_OPS:
            raise _Failure("unknown_op", {"op": op, "handle": self._source_ref(handle) if handle is not None else None})

        if op in _OPAQUE_OPS:
            return self._compile_opaque(handle, row, env, depth)
        if op == "bound_variable":
            return self._compile_bound(handle, row, env)
        if op == "exact_rational":
            return self._compile_rational(handle, row)
        if op == "sum":
            return self._compile_sum(handle, row, env, depth)
        if op == "product":
            return self._compile_product(handle, row, env, depth)
        if op == "negative":
            return self._compile_unary(handle, row, env, depth, "argument")
        if op == "positive_quotient":
            return self._compile_quotient(handle, row, env, depth)
        if op == "analytic_unary":
            return self._compile_analytic_unary(handle, row, env, depth)
        if op == "function_substitution":
            return self._compile_substitution(handle, row, env, depth)
        if op == "definite_integral":
            return self._compile_integral(handle, row, env, depth)
        raise _Failure("unknown_op", {"op": op})

    def _field(self, row, names, required=True):
        value = _first(row, names, _MISSING)
        if value is _MISSING and required:
            raise _Failure("malformed_row", {"missing": list(names)})
        return value

    def _compile_bound(self, handle, row, env):
        name = self._field(row, ("name", "variable", "variable_name"))
        if not isinstance(name, str) or not name:
            raise _Failure("invalid_bound_name")
        if name in env:
            binding = env[name]
            if isinstance(binding, dict) and "__expression__" in binding:
                expected_handle = binding.get("__source_handle__")
                if expected_handle is not None and handle != expected_handle:
                    raise _Failure("same_name_different_handle", {"name": name})
                if expected_handle is None and handle is not None:
                    binding["__source_handle__"] = self._source_ref(handle)
                return _clone(binding["__expression__"])
            # A resolved variable is the already-compiled expression itself.
            # Returning a wrapper here would leave a second, evaluator-specific
            # substitution pass and would make alpha hygiene ambiguous.
            return _clone(binding)
        if self.scope_candidate_handles:
            candidates = self.scope_candidate_handles[-1]
            if name in candidates and candidates[name] != handle:
                raise _Failure("same_name_different_handle", {"name": name})
            candidates[name] = self._source_ref(handle) if handle is not None else None
        result = self._base("bound_variable", handle, row)
        result["name"] = name
        result["resolved"] = False
        result["source_node_reference"] = self._source_ref(handle) if handle is not None else None
        return result

    def _compile_rational(self, handle, row):
        numerator = _first(row, ("numerator", "num"), _MISSING)
        denominator = _first(row, ("denominator", "den"), _MISSING)
        value = _first(row, ("value", "exact_value", "integer"), _MISSING)
        if numerator is _MISSING and value is _MISSING:
            raise _Failure("malformed_exact_rational")
        if numerator is _MISSING:
            numerator = value
        if denominator is _MISSING:
            denominator = 1
        if not isinstance(numerator, int) or isinstance(numerator, bool):
            raise _Failure("invalid_exact_rational_numerator")
        if not isinstance(denominator, int) or isinstance(denominator, bool):
            raise _Failure("invalid_exact_rational_denominator")
        if denominator == 0:
            raise _Failure("zero_exact_rational_denominator")
        result = self._base("exact_rational", handle, row)
        result["numerator"] = _clone(numerator)
        result["denominator"] = _clone(denominator)
        # This is a structural marker only.  It does not claim a computed
        # numeric value for an opaque graph.
        result["is_exact_rational"] = True
        return result

    def _list_field(self, row, names):
        value = self._field(row, names)
        if not isinstance(value, (list, tuple)):
            raise _Failure("malformed_operator_children")
        if not value:
            raise _Failure("empty_operator_children")
        return list(value)

    def _compile_sum(self, handle, row, env, depth):
        children = self._list_field(row, ("arguments", "terms", "children", "operands", "args"))
        result = self._base("sum", handle, row)
        result["terms"] = [self._compile_ref(child, env, depth + 1) for child in children]
        return result

    def _compile_product(self, handle, row, env, depth):
        children = self._list_field(row, ("arguments", "factors", "children", "operands", "args"))
        result = self._base("product", handle, row)
        result["factors"] = [self._compile_ref(child, env, depth + 1) for child in children]
        return result

    def _compile_unary(self, handle, row, env, depth, field_name):
        child = self._field(row, (field_name, "child", "operand", "argument", "arg"))
        result = self._base(_row_op(row), handle, row)
        result["argument"] = self._compile_ref(child, env, depth + 1)
        return result

    def _compile_quotient(self, handle, row, env, depth):
        numerator = self._field(row, ("numerator", "num"))
        denominator = self._field(row, ("denominator", "den"))
        result = self._base("positive_quotient", handle, row)
        result["numerator"] = self._compile_ref(numerator, env, depth + 1)
        result["denominator"] = self._compile_ref(denominator, env, depth + 1)
        return result

    def _compile_analytic_unary(self, handle, row, env, depth):
        child = self._field(row, ("argument", "child", "operand", "arg"))
        result = self._base("analytic_unary", handle, row)
        name = _first(row, ("function", "function_name", "name"), _MISSING)
        if name is _MISSING or not isinstance(name, str) or name.lower() not in _ANALYTIC_UNARY_NAMES:
            raise _Failure("unknown_analytic_unary", {"name": _clone(name) if name is not _MISSING else None})
        result["function"] = _clone(name)
        result["argument"] = self._compile_ref(child, env, depth + 1)
        return result

    def _compile_substitution(self, handle, row, env, depth):
        variable = self._field(row, ("variable", "binder", "parameter"))
        value = self._field(row, ("value", "replacement", "actual"))
        expression = self._field(row, ("expression", "body", "in"))
        variable_name = self._variable_name(variable)
        variable_handle = self._variable_handle(variable)
        compiled_value = self._compile_ref(value, env, depth + 1)
        if len(env) >= _MAX_ENV:
            raise _Failure("environment_limit_exceeded")
        extended = dict(env)
        if variable_name in env:
            previous = env[variable_name]
            if isinstance(previous, dict) and "__source_handle__" in previous:
                previous_handle = previous.get("__source_handle__")
                if previous_handle is not None and variable_handle != previous_handle:
                    raise _Failure("same_name_different_handle", {"name": variable_name})
        extended[variable_name] = {
            "__expression__": compiled_value,
            "__source_handle__": self._source_ref(variable_handle) if variable_handle is not None else None,
        }
        # Function substitution is lexical and source preserving: the binder
        # wrapper itself is represented in metadata while the returned expr is
        # the hygienically compiled body.
        self.scope_candidate_handles.append({})
        try:
            body = self._compile_ref(expression, extended, depth + 1)
        finally:
            self.scope_candidate_handles.pop()
        if isinstance(body, dict):
            body = _clone(body)
            body["substituted_from"] = self._source_ref(handle) if handle is not None else None
            body["substitution_variable"] = variable_name
        return body

    def _variable_name(self, variable):
        if isinstance(variable, dict):
            raise _Failure("unproven_variable_handle")
        # A scalar variable field is a graph handle.  It must resolve to an
        # explicit bound-variable declaration; a same-name string without a
        # node proof is deliberately unresolved rather than treated as a
        # global binder identity.
        if variable is not _MISSING:
            try:
                row = self._lookup(variable)
            except _Failure:
                raise _Failure("unproven_variable_handle")
            if isinstance(row, dict):
                if _row_op(row) != "bound_variable":
                    raise _Failure("invalid_substitution_variable_row")
                name = _first(row, ("name", "variable", "variable_name"), _MISSING)
                if isinstance(name, str) and name:
                    return name
                raise _Failure("invalid_substitution_variable_row")
        raise _Failure("invalid_substitution_variable")

    def _variable_handle(self, variable):
        if isinstance(variable, dict):
            raise _Failure("unproven_variable_handle")
        try:
            row = self._lookup(variable)
        except _Failure:
            raise _Failure("unproven_variable_handle")
        if _row_op(row) != "bound_variable":
            raise _Failure("invalid_substitution_variable_row")
        return variable

    def _compile_integral(self, handle, row, env, depth):
        lower = self._field(row, ("lower", "lower_bound", "a"))
        upper = self._field(row, ("upper", "upper_bound", "b"))
        integrand = self._field(row, ("integrand", "body", "expression"))
        binder = self._field(row, ("variable", "binder", "parameter", "bound_variable"))
        if not isinstance(binder, str):
            raise _Failure("integral_binder_not_string")
        binder_name = binder
        if binder_name != _INTEGRAL_VARIABLE:
            raise _Failure("unexpected_integral_binder", {"name": binder_name})
        compiled_lower = self._compile_ref(lower, env, depth + 1)
        compiled_upper = self._compile_ref(upper, env, depth + 1)
        # The archived entry is one of the unit-integral frontiers.  Keep the
        # endpoint rows structural, but require their exact 0 and 1 shape so a
        # different inlet cannot be mistaken for this Q4 archive.
        if not self._is_exact_integer(compiled_lower, 0) or not self._is_exact_integer(compiled_upper, 1):
            raise _Failure("wrong_integral_endpoints")
        token = self._fresh_binder(handle, binder_name)
        extended = dict(env)
        extended[binder_name] = {
            "__expression__": token,
            "__source_handle__": None,
        }
        self.scope_candidate_handles.append({})
        try:
            compiled_integrand = self._compile_ref(integrand, extended, depth + 1)
        finally:
            self.scope_candidate_handles.pop()
        self.integral_handles.append(self._source_ref(handle) if handle is not None else None)
        result = self._base("definite_integral", handle, row)
        result["lower"] = compiled_lower
        result["upper"] = compiled_upper
        result["binder"] = token
        result["integrand"] = compiled_integrand
        return result

    def _fresh_binder(self, handle, name):
        self.binder_count += 1
        # The token includes a compiler-local ordinal and source handle.  It is
        # never the user-facing source name, so it cannot collide by alpha name.
        ref = self._source_ref(handle) if handle is not None else "inline"
        token = "q4_fresh_binder_" + str(self.binder_count)
        return {
            "kind": "fresh_binder",
            "token": token,
            "source_name": name,
            "source_ref": ref,
        }

    def _compile_opaque(self, handle, row, env, depth):
        op = _row_op(row)
        self._record_atom(handle, row, op)
        if op == "current_original_source_parameter":
            if row.get('source_family') != self.family or row.get('Z_independent') is not True or not isinstance(row.get('name'), str):
                raise _Failure('malformed_source_parameter')
            result = self._base(op, handle, row)
            result["kind"] = "symbolic_source"
            result["symbolic"] = True
            result["numeric_constant"] = False
            return result
        if op == "current_original_unique_analytic_core_component":
            component_family = _first(row, ("source_family", "family"), _MISSING)
            if component_family is _MISSING or not _family_equal(component_family, self.family):
                raise _Failure("component_family_mismatch")
            rho = self._field(row, ("rho",))
            z_variable = self._field(row, ("Z_variable", "z_variable", "Z"))
            provider = _first(row, ("provider", "provider_id", "provider_handle"), _MISSING)
            provider_refs = {
                key: _clone(value)
                for key, value in row.items()
                if isinstance(key, str) and "provider" in key.lower()
            }
            nested_providers = row.get("providers")
            if isinstance(nested_providers, dict):
                for key, value in nested_providers.items():
                    provider_refs.setdefault(str(key), _clone(value))
            if provider is _MISSING and not provider_refs:
                raise _Failure("malformed_analytic_component", {"missing": "provider provenance"})
            if len(provider_refs) < 4:
                raise _Failure("malformed_analytic_component", {"missing": "four provider refs"})
            if provider is _MISSING:
                provider = provider_refs
            quantity = self._field(row, ("quantity",))
            if quantity not in ('Phi', 'Psi', 'mean_Psi'):
                raise _Failure('unknown_component_quantity')
            derivative = _first(
                row,
                (
                    "Z_derivative",
                    "z_derivative",
                    "Z_derivative_order",
                    "z_derivative_order",
                    "ordinary_slow_Z_derivative_order",
                    "ordinary_slow_z_derivative_order",
                    "derivative",
                ),
                _MISSING,
            )
            if type(derivative) is not int or derivative < 0:
                raise _Failure("malformed_analytic_component", {"missing": "Z_derivative"})
            for argument_name, argument_value in (("rho", rho), ("Z_variable", z_variable)):
                if isinstance(argument_value, dict):
                    if _row_op(argument_value) != "bound_variable":
                        raise _Failure("component_argument_not_bound_variable", {"argument": argument_name})
                else:
                    try:
                        argument_row = self._lookup(argument_value)
                    except _Failure:
                        raise _Failure("missing_component_argument", {"argument": argument_name})
                    if _row_op(argument_row) != "bound_variable":
                        raise _Failure("component_argument_not_bound_variable", {"argument": argument_name})
            result = self._base(op, handle, row)
            result["kind"] = "analytic_core_component"
            result["provider"] = _clone(provider)
            result["provider_refs"] = provider_refs
            result["quantity"] = _clone(quantity)
            result["Z_derivative"] = _clone(derivative)
            # Keep the complete archive row and its named provenance fields
            # visible without asking any provider to evaluate anything.
            result["provenance"] = _clone(row)
            for name in (
                "source_family",
                "radial_order",
                "ordinary_slow_Z_derivative_order",
                "ordinary_slow_z_derivative_order",
            ):
                if name in row:
                    result[name] = _clone(row[name])
            result["arguments"] = {
                "rho": self._compile_ref(rho, env, depth + 1),
                "Z_variable": self._compile_ref(z_variable, env, depth + 1),
            }
            result["numeric_solve"] = False
            return result
        raise _Failure("unknown_op", {"op": op})

    def _family_check(self):
        if not _family_is_valid(self.family):
            raise _Failure("family_mismatch", {"source_family": _clone(self.family)})
        if self.family_name and ("p1" in self.family_name or "m15" in self.family_name):
            raise _Failure("family_mismatch", {"source_family": _clone(self.family)})
        archive_family = _first(self.metadata, ("source_family", "family"), _MISSING)
        if archive_family is not _MISSING:
            if not _family_equal(archive_family, self.family):
                raise _Failure("family_mismatch", {"archive_family": _clone(archive_family)})

    def _root(self):
        row = self._node(self.inlet_handle)
        if _row_op(row) != "function_substitution":
            raise _Failure("wrong_outer_operator")
        self.outer_handle = self.inlet_handle
        variable = self._field(row, ("variable", "binder", "parameter"))
        variable_name = self._variable_name(variable)
        if variable_name != _OUTER_VARIABLE:
            raise _Failure("wrong_outer_variable", {"name": variable_name})
        value = self._field(row, ("value", "replacement", "actual"))
        value_row = self._node(value) if not isinstance(value, dict) else value
        if _row_op(value_row) != "exact_rational":
            raise _Failure("outer_value_not_exact_rational")
        value_ir = self._compile_ref(value, {}, 0)
        if not self._is_exact_four(value_ir):
            raise _Failure("outer_value_not_exact_four")
        expression = self._field(row, ("expression", "body", "in"))
        # Compile the value in the outer environment, then extend by the named
        # variable exactly once.  This ordering prevents accidental capture.
        extended = {
            _OUTER_VARIABLE: {
                "__expression__": value_ir,
                "__source_handle__": self._source_ref(variable),
            }
        }
        self.scope_candidate_handles.append({})
        try:
            return self._compile_ref(expression, extended, 0)
        finally:
            self.scope_candidate_handles.pop()

    def _is_exact_four(self, node):
        if not isinstance(node, dict) or node.get("op") != "exact_rational":
            return False
        return node.get("numerator") == 4 and node.get("denominator") == 1

    def _is_exact_integer(self, node, expected):
        if not isinstance(node, dict) or node.get("op") != "exact_rational":
            return False
        return node.get("numerator") == expected and node.get("denominator") == 1

    def compile(self):
        self._family_check()
        ir = self._root()
        # The root must contain the target integral frontier.  An outer
        # substitution that never reaches it is not this archived Q4 entry.
        if not self.integral_handles:
            raise _Failure("missing_target_integral")
        return ir


def _meta_copy(metadata, names):
    result = {}
    for name in names:
        if name in metadata:
            result[name] = _clone(metadata[name])
    return result


def _find_integral(node):
    if isinstance(node, dict):
        if node.get("op") == "definite_integral":
            return node
        for value in node.values():
            found = _find_integral(value)
            if found is not None:
                return found
    elif isinstance(node, list):
        for value in node:
            found = _find_integral(value)
            if found is not None:
                return found
    return None


def _is_exact_four(node):
    return (isinstance(node, dict) and node.get('op') == 'exact_rational'
            and node.get('numerator') == 4 and node.get('denominator') == 1)


def _component_reports(node, current_binder_token=None):
    """Report every core component without interpreting its value."""

    reports = []
    if isinstance(node, dict):
        if node.get("op") == "definite_integral":
            binder = node.get("binder")
            if isinstance(binder, dict):
                current_binder_token = binder.get("token")
        if node.get("kind") == "analytic_core_component":
            rho_argument = node.get("arguments", {}).get("rho")
            is_four_times_fresh = False
            if isinstance(rho_argument, dict) and rho_argument.get("op") == "product":
                factors = rho_argument.get("factors")
                if isinstance(factors, list) and len(factors) == 2:
                    is_four = any(_is_exact_four(item) for item in factors)
                    is_fresh = any(
                        isinstance(item, dict)
                        and item.get("kind") == "fresh_binder"
                        and item.get("token") == current_binder_token
                        for item in factors
                    )
                    is_four_times_fresh = is_four and is_fresh
            reports.append(
                {
                    "source_ref": _clone(node.get("source_ref")),
                    "rho_argument": _clone(rho_argument),
                    "Z_variable_argument": _clone(node.get("arguments", {}).get("Z_variable")),
                    "rho_is_exact_rational_4_times_fresh_binder": is_four_times_fresh,
                    "provider": _clone(node.get("provider")),
                    "provider_refs": _clone(node.get("provider_refs", {})),
                    "quantity": _clone(node.get("quantity")),
                    "Z_derivative": _clone(node.get("Z_derivative")),
                }
            )
        for value in node.values():
            reports.extend(_component_reports(value, current_binder_token))
    elif isinstance(node, list):
        for value in node:
            reports.extend(_component_reports(value, current_binder_token))
    return reports


def _is_integral_with_fresh_binder(node):
    return (
        isinstance(node, dict)
        and node.get("op") == "definite_integral"
        and isinstance(node.get("binder"), dict)
        and node["binder"].get("kind") == "fresh_binder"
    )


def _replacement_report(ir, component_reports=None):
    if component_reports is None:
        component_reports = _component_reports(ir)
    rho_matches = [
        item
        for item in component_reports
        if item.get("rho_is_exact_rational_4_times_fresh_binder") is True
    ]
    integral = _find_integral(ir)
    report = {
        "expression": _clone(ir),
        "is_exact_rational_4_times_integral_fresh_binder": bool(rho_matches),
        "is_component_rho_exact_rational_4_times_integral_fresh_binder": bool(rho_matches),
        "exact_rational_4_source_ref": None,
        "integral_source_ref": _clone(integral.get("source_ref")) if integral is not None else None,
    }
    if rho_matches:
        factors = rho_matches[0].get("rho_argument", {}).get("factors", [])
        for factor in factors if isinstance(factors, list) else []:
            if _is_exact_four(factor):
                report["exact_rational_4_source_ref"] = _clone(factor.get("source_ref"))
                break
    return report


def _unresolved_result(compiler, family, reason, detail=None):
    metadata = compiler.metadata if compiler is not None else {}
    source_family = _clone(family)
    reasons = []
    if reason:
        reasons.append({"code": reason, "detail": _clone(detail) if detail is not None else None})
    if compiler is not None:
        for item in compiler.reasons:
            reasons.append(_clone(item))
    source_atoms = _clone(compiler.source_atoms) if compiler is not None else []
    source_map = _clone(compiler.source_atom_map) if compiler is not None else {}
    result = {
        "status": "unresolved",
        "source_family": source_family,
        "ir": None,
        "source_atoms": source_atoms,
        "source_atoms_mapping": source_map,
        "report": {
            "source_atoms": source_atoms,
            "source_atoms_mapping": source_map,
            "inner_replacement": None,
            "core_components": [],
        },
        "errors": reasons,
        "scope": {
            "capture_safe_for_this_IR": False,
            "graph_runtime_semantics_proved": False,
            "global_string_to_node_identity_proved": False,
            "numeric_field_value_computed": False,
            "old_errors_removed": False,
            "terminal_closure_proved": False,
        },
    }
    result.update(_meta_copy(metadata, ("bounds", "old_errors", "numeric", "numeric_field", "numeric_field_value")))
    if "errors" in metadata:
        result["archive_errors"] = _clone(metadata["errors"])
    return result


def build(graph, inlet_handle, family):
    """Lower one archived Q4 inlet into a source-preserving hygienic IR.

    ``graph`` may be ``{"exact_graph_nodes": [row, ...], ...metadata...}``
    (the archived shape), ``{"nodes": {handle: row}, ...metadata...}``, a
    graph with a list of rows carrying ``handle``/``id``, or a direct
    handle-to-row mapping.  The result is always JSON-compatible.  Any unknown, malformed,
    cyclic, over-large, wrong-family, or wrong-entry graph returns
    ``status == "unresolved"`` with no partial claim of runtime semantics.
    """

    compiler = None
    try:
        compiler = _Compiler(graph, inlet_handle, family)
        ir = compiler.compile()
        source_atoms = _clone(compiler.source_atoms)
        source_map = _clone(compiler.source_atom_map)
        component_reports = _component_reports(ir)
        result = {
            "status": "lowered",
            "source_family": _clone(family),
            "ir": _clone(ir),
            "source_atoms": source_atoms,
            "source_atoms_mapping": source_map,
            "report": {
                "source_atoms": source_atoms,
                "source_atoms_mapping": source_map,
                "inner_replacement": _replacement_report(ir, component_reports),
                "core_components": component_reports,
            },
            "scope": {
                "capture_safe_for_this_IR": True,
                "graph_runtime_semantics_proved": False,
                "global_string_to_node_identity_proved": False,
                "numeric_field_value_computed": False,
                "old_errors_removed": False,
                "terminal_closure_proved": False,
            },
        }
        result.update(_meta_copy(compiler.metadata, ("bounds", "old_errors", "numeric", "numeric_field", "numeric_field_value")))
        if "errors" in compiler.metadata:
            result["archive_errors"] = _clone(compiler.metadata["errors"])
        return result
    except _Failure as exc:
        return _unresolved_result(compiler, family, exc.reason, exc.detail)
    except (KeyError, TypeError, ValueError, RecursionError) as exc:
        # A malformed archive must never escape as an evaluator exception.
        return _unresolved_result(compiler, family, "invalid_graph_value", str(exc))


__all__ = ["build"]
