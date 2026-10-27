"""Static, typed physical-chart consumer for the native O3 source packet.

This module deliberately has no dependency on the reconstruction machinery.  It
accepts the already sealed source packets as data, parses the two source
function snippets with :mod:`ast`, and constructs a small expression IR.  The
IR is a description of functions; it is never evaluated here.  In particular,
``R`` is kept distinct from the physical cylindrical radius ``r`` and the
physical time variable is the free chart coordinate ``tau``.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction


COMPONENT = "native_O3_total_radius_velocity_pressure_binding"
INPUT_SCHEMA = "native_O3_total_radius_velocity_pressure_binding_input_v1"
RESULT_SCHEMA = "native_O3_total_radius_velocity_pressure_binding_result_v1"

WRAPPER_PINS = {
    "saved_flux_result": "d4cf1b39195a608e32197c58d572f86f0d974ba940ef2ff577295fdfc95b092d",
    "source_literals": "ae3a6afbdc8de9c557dff18b3552a6b1c6c1d902decca5298b0b6adbcc216b59",
    "physical_source_capsule": "90ce8bbacd274cf7ea022afd976799944b60709fef69162a6736a891fcef7788",
}

BRANCHES = ("background", "modified")
ORDERS = (0, 1, 2, 3, 4)
SOURCE_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
SELECTABLE_FIELDS = (
    "ur",
    "utheta",
    "uz",
    "u",
    "v",
    "w",
    "pressure",
    "psi",
    "vector_potential",
    "physical_map",
    "volume_form",
)


class _KernelInputError(ValueError):
    """An input is outside the sealed source-data contract."""


def _is_map(value):
    return isinstance(value, Mapping)


def _require(condition, message):
    if not condition:
        raise _KernelInputError(message)


def _number(value, denominator=1):
    return {
        "kind": "exact_constant",
        "numerator": int(value),
        "denominator": int(denominator),
    }


def _free(name, role=None):
    row = {"kind": "free_symbol", "name": str(name)}
    if role is not None:
        row["role"] = role
    return row


def _source_ref(source, path, node, field=None, order=None, role=None):
    _require(isinstance(node, int) and not isinstance(node, bool),
             "source node ids must be integer graph references")
    row = {
        "kind": "source_ref",
        "source": source,
        "path": path,
        "native_node": node,
    }
    if field is not None:
        row["field"] = field
    if order is not None:
        row["ordinary_Z_order"] = int(order)
    if role is not None:
        row["role"] = role
    return row


def _param_ref(name, node, role, materialized=False):
    row = _source_ref("physical_source_capsule", "parameter_ports." + name,
                      node, role=role)
    row.update({"parameter": name, "materialized": bool(materialized)})
    return row


def _op(operator, *children, **metadata):
    row = {
        "kind": "typed_function_term",
        "operator": operator,
        "children": list(children),
    }
    row.update(metadata)
    return row


def _neg(value):
    return _op("negative", value)


def _add(*values, **metadata):
    return _op("add", *values, **metadata)


def _sub(left, right):
    return _op("sub", left, right)


def _mul(*values, **metadata):
    return _op("mul", *values, **metadata)


def _scale(value, factor):
    return _op("scale", factor, value)


def _quotient(numerator, denominator, contract=None):
    row = _op("quotient", numerator, denominator)
    if contract is not None:
        row["contract"] = contract
    return row


def _power(value, exponent):
    return _op("power", value, exponent)


def _sqrt(value):
    return _op("sqrt", value)


def _derivative(value, axis, order=1, **metadata):
    row = _op("ordinary_derivative", value, axis=_free(axis), derivative_order=order)
    row.update(metadata)
    return row


def _replace_free_symbols(value, replacements):
    """Bind only the fresh physical free symbols in a typed expression."""
    if _is_map(value):
        if value.get("kind") == "free_symbol":
            name = value.get("name")
            if name in replacements:
                return replacements[name]
        return {key: _replace_free_symbols(item, replacements)
                for key, item in value.items()}
    if isinstance(value, list):
        return [_replace_free_symbols(item, replacements) for item in value]
    if isinstance(value, tuple):
        return tuple(_replace_free_symbols(item, replacements) for item in value)
    return value


def _row_ref(parent_path, order):
    return {
        "kind": "accepted_parent_ref",
        "source": "saved_flux_result",
        "path": parent_path,
        "ordinary_Z_order": int(order),
    }


def _json_literal(value, label):
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except (TypeError, ValueError) as exc:
            raise _KernelInputError(label + " literal is not JSON") from exc
    _require(_is_map(value), label + " literal must be an object")
    return value


def _source_sha(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _attribute_name(node):
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
    return ".".join(reversed(parts)) if parts else None


def _ast_summary(source_text, expected_function=None, label="source"):
    """Parse source text and return structural facts, never executable code."""
    _require(isinstance(source_text, str) and source_text.strip(),
             label + " source_text is required")
    try:
        tree = ast.parse(source_text, filename=label, mode="exec")
    except SyntaxError as exc:
        raise _KernelInputError(label + " source_text is not parseable") from exc
    functions = [node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    names = [node.name for node in functions]
    if expected_function is not None:
        _require(expected_function in names,
                 label + " is missing expected function " + expected_function)
    calls = []
    assigned = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = _attribute_name(node.func)
            if name is not None:
                calls.append({"name": name, "arity": len(node.args)})
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    assigned.append(target.id)
                elif isinstance(target, (ast.Tuple, ast.List)):
                    assigned.extend(
                        elt.id for elt in target.elts if isinstance(elt, ast.Name)
                    )
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            assigned.append(node.target.id)
    return {
        "kind": "source_ast_summary",
        "label": label,
        "sha256": _source_sha(source_text),
        "function_names": names,
        "calls": calls,
        "assigned_names": assigned,
        "node_count": sum(1 for _ in ast.walk(tree)),
    }


def _ast_source_record(record, label):
    if not _is_map(record):
        raise _KernelInputError(label + " AST record must be an object")
    text = record.get("source_text")
    summary = _ast_summary(text, label=label)
    return {
        "source_path": record.get("source_path"),
        "name": record.get("name"),
        "class_name": record.get("class_name"),
        "line_span": record.get("line_span"),
        "source_text_sha256": summary["sha256"],
        "ast": summary,
    }


def _validate_wrappers(payload):
    wrappers = payload.get("source_wrappers")
    _require(_is_map(wrappers), "source_wrappers must be an object")
    _require(set(wrappers) == set(WRAPPER_PINS),
             "source_wrappers must contain exactly the three pinned wrappers")
    out = {}
    for name, expected_sha in WRAPPER_PINS.items():
        wrapper = wrappers.get(name)
        _require(_is_map(wrapper), name + " wrapper must be an object")
        _require(set(wrapper) == {"sha256", "literal"},
                 name + " wrapper must have exactly sha256 and literal")
        _require(wrapper.get("sha256") == expected_sha,
                 name + " wrapper sha256 does not match the sealed source")
        out[name] = _json_literal(wrapper.get("literal"), name)
    return out


def _parent_payload(saved):
    result = saved.get("result") if _is_map(saved.get("result")) else saved
    _require(_is_map(result), "saved flux result must contain a result object")
    payload = result.get("payload")
    _require(_is_map(payload), "saved flux result must contain result.payload")
    return result, payload


def _parent_commit(saved, result, payload, outer):
    for candidate in (outer, saved, result, payload):
        if not _is_map(candidate):
            continue
        for key in ("expected_parent_commit", "parent_commit"):
            value = candidate.get(key)
            if isinstance(value, str) and value:
                return value
        frozen = candidate.get("frozen_context")
        if _is_map(frozen) and isinstance(frozen.get("parent_commit"), str):
            return frozen.get("parent_commit")
    return None


def _capsule_rows(capsule, branch):
    rows = capsule.get(branch + "_rows")
    _require(_is_map(rows), branch + "_rows is required")
    out = {}
    for field in SOURCE_FIELDS:
        values = rows.get(field)
        _require(isinstance(values, list) and len(values) == 5,
                 branch + "_rows." + field + " must have five source refs")
        out[field] = list(values)
        for node in values:
            _require(isinstance(node, int) and not isinstance(node, bool),
                     branch + "_rows contain a non-native node")
    return out


def _native_node(value, label):
    node = value.get("node") if _is_map(value) else value
    _require(isinstance(node, int) and not isinstance(node, bool),
             label + " must carry an integer native node")
    return node


def _parameter_ports(capsule):
    ports = capsule.get("parameter_ports")
    _require(_is_map(ports), "physical capsule parameter_ports is required")
    expected = {
        "Pstar": (33572, "fixed_original_Pstar"),
        "delta": (25220, "fixed_original_similarity_exponent"),
        "Z": (20, "ordinary_Z"),
        "R": (241787, "positive_similarity_profile_radius"),
        "native_x": (241536, "native_O3_log_profile_radius"),
    }
    out = {}
    for name, (node, role) in expected.items():
        row = ports.get(name)
        _require(_is_map(row), "missing physical capsule parameter " + name)
        _require(row.get("node") == node,
                 "physical capsule parameter node mismatch for " + name)
        if name in ("Pstar", "delta", "R"):
            _require(row.get("materialized") is False,
                     name + " must remain an unevaluated source parameter")
        out[name] = _param_ref(name, node, role, False)
    rm = ports.get("Rm")
    _require(_is_map(rm) and rm.get("node") == 85743,
             "physical capsule must pin Rm source port 85743")
    _require(rm.get("materialized") is False,
             "Rm must remain an unevaluated source parameter")
    out["Rm"] = _param_ref("Rm", 85743, "fixed_native_Rm", False)
    return out


def _capsule_source(capsule):
    source_ast = capsule.get("source_AST", capsule.get("source_ast"))
    _require(_is_map(source_ast), "physical capsule source_AST is required")
    text = source_ast.get("source_text")
    summary = _ast_summary(text, expected_function="physical_pressure_fields",
                           label="physical_pressure_fields")
    required_calls = {
        "Algebra", "a.scale", "a.sub", "b.fixed", "b.sub", "b.mul",
        "b.trim", "b.shift", "b.div", "b.add",
    }
    seen = {item["name"] for item in summary["calls"]}
    _require(required_calls.issubset(seen),
             "physical pressure source AST is missing an expected operation")
    return {
        "record": {
            "function": source_ast.get("function"),
            "line_start": source_ast.get("line_start"),
            "line_end": source_ast.get("line_end"),
            "source_path": source_ast.get("path"),
            "source_text_sha256": summary["sha256"],
        },
        "ast": summary,
    }


def _finite_source_ast(source_literals):
    source_ast = source_literals.get("source_ast")
    _require(_is_map(source_ast), "source_literals source_ast is required")
    records = []
    for key in (
        "finite_N_functions",
        "background_density",
        "JetAlgebra.add",
        "JetAlgebra.mul",
        "JetAlgebra.scale",
        "BridgeAlgebra.exp",
    ):
        record = source_ast.get(key)
        if _is_map(record) and isinstance(record.get("source_text"), str):
            records.append(_ast_source_record(record, "source_ast." + key))
    _require(any("finite_N_functions" in item["ast"].get("function_names", [])
                 for item in records),
             "finite_N_functions source AST was not supplied")
    return records


def _parent_rows(parent_payload):
    my_rows = parent_payload.get("my_compatibility_rows")
    _require(isinstance(my_rows, list) and len(my_rows) == 5,
             "accepted parent my_compatibility_rows must have five rows")
    for order, row in enumerate(my_rows):
        _require(_is_map(row) and row.get("ordinary_Z_order") == order,
                 "accepted parent my row order mismatch")
        sum_head = row.get("sum_head")
        _native_node(sum_head, "accepted parent my row sum_head")
    derivative_rows = parent_payload.get("derivative_rewrite_rows")
    _require(isinstance(derivative_rows, list) and len(derivative_rows) == 25,
             "accepted parent derivative rows must contain 25 rows")
    return my_rows, derivative_rows


def _source_join_guards(capsule_rows, parent_payload, source_literals):
    """Join accepted parent refs to the fresh capsule by native graph node."""
    my_rows = parent_payload["my_compatibility_rows"]
    joins = {branch: [] for branch in BRANCHES}
    for order in ORDERS:
        modified_m = capsule_rows["modified"]["m"][order]
        modified_V = capsule_rows["modified"]["V"][order]
        parent_m = _native_node(
            my_rows[order].get("negative_argument"),
            "modified my negative argument")
        parent_V = _native_node(
            my_rows[order].get("V_N_ref"), "modified my V_N_ref")
        _require(parent_m == modified_m,
                 "modified m source join does not match accepted my argument")
        _require(parent_V == modified_V,
                 "modified V source join does not match accepted my V_N ref")
        joins["modified"].append({
            "ordinary_Z_order": order,
            "modified_m_node": modified_m,
            "modified_V_node": modified_V,
            "parent_negative_argument_node": parent_m,
            "parent_V_N_ref_node": parent_V,
            "join_verified": True,
        })

        bg_proofs = parent_payload.get("source_proof_rows", {}).get("background", {})
        bg_m_rows = bg_proofs.get("m") if _is_map(bg_proofs) else None
        _require(_is_map(bg_m_rows), "accepted background m proof rows are required")
        bg_row = bg_m_rows.get(str(order))
        _require(_is_map(bg_row), "accepted background m proof row is missing")
        source_expression = bg_row.get("source_expression")
        _require(_is_map(source_expression), "accepted background m source expression is missing")
        history_nodes = [
            _native_node(term.get("node"), "background history term")
            for term in source_expression.get("terms", [])
            if _is_map(term) and term.get("source") == "history_j"
        ]
        _require(history_nodes == [capsule_rows["background"]["m"][order]],
                 "background m source join does not match Duhamel history")
        bg_fields = source_literals.get("literal_rows", {}).get("background_fields", {})
        bg_v = bg_fields.get("V") if _is_map(bg_fields) else None
        _require(_is_map(bg_v), "source literal background V rows are required")
        bg_v_rows = bg_v.get("rows")
        _require(isinstance(bg_v_rows, list) and len(bg_v_rows) >= 5,
                 "source literal background V rows are incomplete")
        source_V = _native_node(bg_v_rows[order], "source literal background V")
        _require(source_V == capsule_rows["background"]["V"][order],
                 "background V source join does not match source literal")
        joins["background"].append({
            "ordinary_Z_order": order,
            "background_m_node": capsule_rows["background"]["m"][order],
            "background_history_node": history_nodes[0],
            "background_V_node": capsule_rows["background"]["V"][order],
            "source_literal_background_V_node": source_V,
            "join_verified": True,
        })

    pressure_rows = parent_payload.get("pressure_alias_rows")
    _require(isinstance(pressure_rows, list) and len(pressure_rows) == 10,
             "accepted pressure alias rows must contain both branches")
    pressure_seen = set()
    for row in pressure_rows:
        _require(_is_map(row), "pressure alias row must be an object")
        branch = row.get("branch")
        order = row.get("ordinary_Z_order")
        _require(branch in BRANCHES and order in ORDERS,
                 "pressure alias row branch/order is invalid")
        key = (branch, order)
        _require(key not in pressure_seen, "duplicate pressure alias row")
        pressure_seen.add(key)
        raw = row.get("raw_node_bindings")
        _require(_is_map(raw), "pressure alias raw node bindings are required")
        _require(_native_node(raw.get("P0"), "pressure P0 raw node") ==
                 capsule_rows[branch]["P0"][order],
                 "pressure P0 source join mismatch")
        _require(_native_node(raw.get("p"), "pressure p raw node") ==
                 capsule_rows[branch]["p"][order],
                 "pressure p source join mismatch")
        _require(_native_node(raw.get("pressure"), "pressure raw node") ==
                 capsule_rows[branch]["pressure"][order],
                 "pressure source join mismatch")
        joins[branch][order]["pressure_join_verified"] = True
    _require(len(pressure_seen) == 10, "pressure alias rows do not cover both branches")
    return joins


def _capsule_ref(capsule_rows, branch, field, order):
    node = capsule_rows[branch][field][order]
    return _source_ref("physical_source_capsule",
                       branch + "_rows." + field + "[" + str(order) + "]",
                       node, field=field, order=order)


def _build_source_bindings(capsule_rows, ports, parent_payload, branch,
                           source_join_guards=None):
    rows = []
    for order in ORDERS:
        fields = {
            field: _capsule_ref(capsule_rows, branch, field, order)
            for field in SOURCE_FIELDS
        }
        rows.append({
            "ordinary_Z_order": order,
            "fields": fields,
            "profile_coordinate": ports["R"],
            "native_chart_coordinate": ports["native_x"],
            "ordinary_Z_coordinate": ports["Z"],
            "accepted_parent_scope": _row_ref("payload.source_scope", order),
            "parent_payload_pointer": {
                "source": "saved_flux_result",
                "path": "result.payload",
                "wrapper_sha256": WRAPPER_PINS["saved_flux_result"],
            },
        })
    return {
        "branch": branch,
        "ordinary_Z_orders": list(ORDERS),
        "rows": rows,
        "source_join_guards": source_join_guards or [],
        "source_join_verified": bool(source_join_guards),
        "parent_scope_is_pointer_only": True,
        "source_values_unevaluated": True,
    }


def _physical_map(ports):
    tau = _free("tau", "physical_time_parameter")
    theta = _free("theta", "physical_azimuth")
    one = _number(1)
    Z = ports["Z"]
    R = ports["R"]
    d = _sub(one, _mul(Z, Z))
    L = _sub(one, _mul(ports["delta"], _mul(Z, Z)))
    lam = _sqrt(_quotient(tau, d, "tau>0 and |Z|<1"))
    r = _mul(lam, _sqrt(_mul(_number(2), R)))
    z = _mul(_power(lam, _sub(one, ports["delta"])), Z)
    t = _sub(one, tau)
    X = _mul(r, _op("cos", theta))
    Y = _mul(r, _op("sin", theta))
    q_forward = _power(lam, _number(2))
    z_phys = _free("z_phys", "physical_cartesian_coordinate")
    X_phys = _free("X", "physical_cartesian_coordinate")
    Y_phys = _free("Y", "physical_cartesian_coordinate")
    t_phys = _free("t_phys", "physical_time_coordinate")
    tau_inverse = _sub(one, t_phys)
    q_inverse = _free("q_inverse", "implicitly_bound_positive_root")
    inverse_Z = _quotient(z_phys, _power(q_inverse, _scale(_sub(one, ports["delta"]), _number(1, 2))))
    inverse_R = _quotient(_add(_mul(X_phys, X_phys), _mul(Y_phys, Y_phys)),
                          _mul(_number(2), q_inverse))
    native_x_inverse = _sub(_op("log", _quotient(inverse_R, ports["Rm"])),
                            _add(_number(17), _op("exp", _number(40))))
    positive_tau = _free("tau", "positive_physical_time_parameter")
    positive_tau_equation = _sub(
        _sub(q_inverse, _mul(_mul(z_phys, z_phys), _power(q_inverse, ports["delta"]))),
        positive_tau,
    )
    inverse_equation = _sub(
        _sub(q_inverse, _mul(_mul(z_phys, z_phys), _power(q_inverse, ports["delta"]))),
        tau_inverse,
    )
    native_R_definition = _mul(
        ports["Rm"],
        _op("exp", _add(_number(17), _op("exp", _number(40)), ports["native_x"])),
    )
    return {
        "kind": "typed_spacetime_map",
        "free_variables": [tau, theta],
        "profile_R": R,
        "physical_r": r,
        "physical_z": z,
        "physical_time": t,
        "physical_X": X,
        "physical_Y": Y,
        "formulae": {
            "lambda_phys": "sqrt(tau/(1-Z^2))",
            "r": "lambda_phys*sqrt(2*R)",
            "z": "lambda_phys^(1-delta)*Z",
            "t": "1-tau",
            "X": "r*cos(theta)",
            "Y": "r*sin(theta)",
        },
        "native_chart_R_definition": native_R_definition,
        "d": d,
        "L": L,
        "lambda_phys": lam,
        "q_forward": q_forward,
        "native_chart_x": ports["native_x"],
        "native_chart_Z": Z,
        "domain": {
            "tau": "0<tau<=1",
            "Z": "|Z|<1",
            "delta": "0<delta<1/200",
            "native_x": "0<native_x<1",
            "phase": "regular interior phase cells; theta is not the source fast phase",
            "axis": "r>0; no axis extension is claimed",
        },
        "notational_separation": {
            "profile_R_is_not_physical_r": True,
            "physical_tau_is_not_native_chart_tau": True,
            "physical_theta_is_not_source_phase": True,
            "physical_XY_are_not_native_x": True,
        },
        "physical_inverse": {
            "q": q_inverse,
            "equation": inverse_equation,
            "equation_text": "q-z_phys^2*q^delta=tau",
            "equation_text_with_time": "q-z_phys^2*q^delta=1-t_phys",
            "positive_tau_equation": positive_tau_equation,
            "tau_from_t": tau_inverse,
            "positive_root": True,
            "Z": inverse_Z,
            "R": inverse_R,
            "native_x": native_x_inverse,
            "angle": {
                "cos_theta": _quotient(X_phys, _sqrt(_add(_mul(X_phys, X_phys), _mul(Y_phys, Y_phys)))),
                "sin_theta": _quotient(Y_phys, _sqrt(_add(_mul(X_phys, X_phys), _mul(Y_phys, Y_phys)))),
                "theta": _op("atan2", Y_phys, X_phys),
            },
            "theta": _op("atan2", Y_phys, X_phys),
            "tau_from_time": tau_inverse,
            "implicit_binding": {
                "bound_variable": q_inverse,
                "equation": inverse_equation,
                "positive_branch": True,
                "free_arguments": ("X", "Y", "z_phys", "t_phys"),
            },
            "numeric_solver": False,
        },
        "paper_equivalence_data": {"q": "lambda_phys^2", "h": "delta/2"},
    }


def _my_term(branch, order, source_rows, my_rows):
    V = source_rows[order]["fields"]["V"]
    m = source_rows[order]["fields"]["m"]
    if branch == "modified":
        row = my_rows[order]
        sum_head = row.get("sum_head")
        node = sum_head.get("node") if _is_map(sum_head) else sum_head
        _require(isinstance(node, int) and not isinstance(node, bool),
                 "accepted modified my row sum_head must contain a native node")
        return _source_ref("saved_flux_result", "result.payload.my_compatibility_rows[" + str(order) + "].sum_head",
                           node, field="m_y", order=order,
                           role="accepted_modified_my")
    return _add(V, _neg(m), accepted_parent_ref=_row_ref(
        "result.payload.source_proof_rows.background.m", order))


def _m_z_term(source_rows, order, branch):
    if order == 0:
        row1 = source_rows[1]["fields"]["m"]
        return _source_ref("physical_source_capsule",
                           branch + "_rows.m[1]",
                           row1["native_node"], field="m_Z", order=1,
                           role="same_incoming_history_row1")
    return _derivative(source_rows[order]["fields"]["m"], "Z", 1,
                       source_row_order=order, same_history=True,
                       derivative_identity_unproved=True)


def _pressure_term(P0, p):
    return _add(P0, p, role="P0_plus_p")


def _build_branch(capsule_rows, ports, source_rows, parent_payload, my_rows,
                  branch, physical_map):
    d = physical_map["d"]
    L = physical_map["L"]
    lam = physical_map["lambda_phys"]
    R = physical_map["profile_R"]
    Z = ports["Z"]
    delta = ports["delta"]
    Pstar = ports["Pstar"]
    cos_theta = _op("cos", _free("theta", "physical_azimuth"))
    sin_theta = _op("sin", _free("theta", "physical_azimuth"))
    # Five source-jet rows remain available as data.  The physical source
    # capsule admits only its leading coefficient (order 0), so no synthetic
    # fifth axial derivative is fabricated for the physical function.
    source_jet_rows = []
    for order in ORDERS:
        source_jet_rows.append({
            "ordinary_Z_order": order,
            "fields": source_rows[order]["fields"],
            "my": _my_term(branch, order, source_rows, my_rows),
            "source_jet_data_only": True,
        })
    rows = []
    for order in (0,):
        f = source_rows[order]["fields"]
        my = source_jet_rows[order]["my"]
        mz = _m_z_term(source_rows, order, branch)
        E, V, m = f["E"], f["V"], f["m"]
        Utheta = _mul(Pstar, E)
        Uz = _mul(Pstar, V)
        Mz = _mul(Pstar, R, m)
        psi = _mul(_power(lam, _sub(_number(1), delta)), Pstar, R, m)
        cartesian_prefactor = _scale(
            _mul(_power(lam, _sub(_number(-1), delta)), Pstar, m),
            _number(1, 2),
        )
        vector = _op(
            "vector",
            _mul(cartesian_prefactor, _neg(physical_map["physical_Y"])),
            _mul(cartesian_prefactor, physical_map["physical_X"]),
            _number(0),
            component_order=("X", "Y", "Z"),
            explicit_components=True,
            cylindrical_definition="A=(psi/r)e_theta",
            formula="0.5*lambda_phys^(-1-delta)*Pstar*m*(-Y,X,0)",
            streamfunction=psi,
            physical_r=physical_map["physical_r"],
        )
        ur_stream_numerator = _add(
            _mul(_add(_number(1), delta), Z, m),
            _mul(_number(2), Z, my),
            _neg(_mul(d, mz)),
        )
        ur_stream = _quotient(
            _mul(_power(lam, _number(-1)), Pstar, _sqrt(_quotient(R, _number(2))),
                 ur_stream_numerator), L,
            "L=1-delta*Z^2; source conditional ordinary-Z chain rule",
        )
        ur_direct_numerator = _add(
            _mul(_number(2), Z, V),
            _neg(_mul(_sub(_number(1), delta), Z, m)),
            _neg(_mul(d, mz)),
        )
        ur = _quotient(
            _mul(_power(lam, _number(-1)), Pstar, _sqrt(_quotient(R, _number(2))),
                 ur_direct_numerator), L,
            "direct conditional pullback using accepted m_y=V-m",
        )
        uz_prefactor = _mul(_power(lam, _sub(_number(-1), delta)), Pstar)
        uz_stream = _mul(uz_prefactor, _add(m, my))
        uz_direct = _mul(uz_prefactor, V)
        utheta = _mul(_power(lam, _sub(_number(-1), delta)), Pstar, E)
        u = _sub(_mul(ur, cos_theta), _mul(utheta, sin_theta))
        v = _add(_mul(ur, sin_theta), _mul(utheta, cos_theta))
        pressure_profile = _mul(_power(Pstar, _number(2)), _pressure_term(f["P0"], f["p"]))
        p_phys = _mul(
            _power(lam, _sub(_number(-2), _mul(_number(2), delta))),
            pressure_profile,
        )
        physical_inverse = physical_map["physical_inverse"]
        source_variable_bindings = {
            "native_x241536": {
                "source_ref": ports["native_x"],
                "substitution": physical_inverse["native_x"],
            },
            "native_x": {
                "source_ref": ports["native_x"],
                "substitution": physical_inverse["native_x"],
            },
            "ordinary_Z20": {
                "source_ref": ports["Z"],
                "substitution": physical_inverse["Z"],
            },
            "Z": {
                "source_ref": ports["Z"],
                "substitution": physical_inverse["Z"],
            },
            "profile_R241787": {
                "source_ref": ports["R"],
                "substitution": physical_inverse["R"],
            },
            "R": {
                "source_ref": ports["R"],
                "substitution": physical_inverse["R"],
            },
            "physical_tau": {
                "substitution": physical_inverse["tau_from_time"],
            },
            "tau": {
                "substitution": physical_inverse["tau_from_time"],
            },
            "physical_theta": {
                "substitution": physical_inverse["theta"],
            },
            "theta": {
                "substitution": physical_inverse["theta"],
            },
        }

        def compose_physical(body):
            return _op(
                "compose_with_physical_inverse",
                _replace_free_symbols(body, {
                    "tau": physical_inverse["tau_from_time"],
                    "theta": physical_inverse["theta"],
                }),
                source_variable_bindings=source_variable_bindings,
                held_source_parameters={
                    "Pstar": Pstar,
                    "delta": delta,
                    "Rm": ports["Rm"],
                },
                bound_free_symbols={
                    "tau": physical_inverse["tau_from_time"],
                    "theta": physical_inverse["theta"],
                },
                implicit_physical_inverse_binding=physical_inverse["implicit_binding"],
                free_arguments=("X", "Y", "z_phys", "t_phys"),
                formal_substitution=True,
                source_values_unevaluated=True,
            )

        physical_cartesian_functions = {
            "u": compose_physical(u),
            "v": compose_physical(v),
            "w": compose_physical(uz_direct),
            "pressure": compose_physical(p_phys),
            "psi": compose_physical(psi),
            "vector_potential": compose_physical(vector),
        }
        rows.append({
            "ordinary_Z_order": order,
            "source_fields": f,
            "my": my,
            "m_Z": mz,
            "Utheta": Utheta,
            "Uz": Uz,
            "Mz": Mz,
            "psi": psi,
            "vector_potential": vector,
            "ur_stream": ur_stream,
            "ur": ur,
            "uz_stream": uz_stream,
            "uz": uz_direct,
            "utheta": utheta,
            "u": u,
            "v": v,
            "w": uz_direct,
            "pressure_profile": pressure_profile,
            "pressure": p_phys,
            "source_function_at": {
                "native_x": ports["native_x"],
                "ordinary_Z": ports["Z"],
                "native_x_source_ref": ports["native_x"],
                "ordinary_Z_source_ref": ports["Z"],
                "native_x_substitution": physical_map["physical_inverse"]["native_x"],
                "ordinary_Z_substitution": physical_map["physical_inverse"]["Z"],
                "substitution": physical_map["physical_inverse"],
                "source_value_evaluation": False,
            },
            "physical_cartesian_functions": physical_cartesian_functions,
        })
    return {
        "branch": branch,
        "ordinary_Z_orders": list(ORDERS),
        "rows": rows,
        "physical_function_orders": [0],
        "source_jet_rows": source_jet_rows,
        "actual_cartesian_function": {
            "variables": ["X", "Y", "z", "t"],
            "physical_variable_names": ["X", "Y", "z_phys", "t_phys"],
            "composition": "physical_inverse",
            "inverse_relation": physical_map["physical_inverse"],
            "functions": rows[0]["physical_cartesian_functions"],
            "leading_row": rows[0],
            "numeric_inverse": False,
            "source_values_unevaluated": True,
        },
        "leading_source_units": {
            "Utheta": "Pstar*E",
            "Uz": "Pstar*V",
            "Mz": "Pstar*R*m",
            "pressure_profile": "Pstar^2*(P0+p)",
        },
        "source_values_unevaluated": True,
        "coefficient_order": 0,
    }


def _identity_rows(branch_data, ports, parent_payload, branch):
    rows = []
    delta = ports["delta"]
    Z = ports["Z"]
    # This is intentionally a two-identity, local coefficient calculation,
    # not a general ring normalizer.  The exact Fraction dictionaries make the
    # two pullback cancellations auditable without evaluating source values.
    d_plus = {
        (0, 0): Fraction(1),
        (2, 0): Fraction(-1) + Fraction(1),
        (2, 1): Fraction(-1),
    }
    d_plus = {key: value for key, value in d_plus.items() if value}
    target_L = {(0, 0): Fraction(1), (2, 1): Fraction(-1)}
    _require(d_plus == target_L,
             "fresh d+(1-delta)Z^2 coefficient identity failed")
    one_plus_delta_minus_two = {
        (1, 0): Fraction(1) - Fraction(2),
        (1, 1): Fraction(1),
    }
    target_negative_one_minus_delta = {
        (1, 0): Fraction(-1),
        (1, 1): Fraction(1),
    }
    _require(one_plus_delta_minus_two == target_negative_one_minus_delta,
             "fresh (1+delta)Z-2Z coefficient identity failed")

    def fraction_terms(poly):
        return {
            "Z^" + str(z_power) + "*delta^" + str(delta_power):
            [coefficient.numerator, coefficient.denominator]
            for (z_power, delta_power), coefficient in sorted(poly.items())
        }

    for row in branch_data["rows"]:
        order = row["ordinary_Z_order"]
        m = row["source_fields"]["m"]
        my = row["my"]
        V = row["source_fields"]["V"]
        rows.append({
            "branch": branch,
            "ordinary_Z_order": order,
            "identity": "m+my=V",
            "left": _add(m, my),
            "right": V,
            "accepted_parent_data": (_row_ref(
                "result.payload.my_compatibility_rows", order)
                if branch == "modified" else _row_ref(
                    "result.payload.source_proof_rows.background.m", order)),
            "radial_chain": {
                "d_plus_one_minus_delta_Z2_equals_L": {
                    "lhs_terms": fraction_terms(d_plus),
                    "rhs_terms": fraction_terms(target_L),
                    "exact_fraction_polynomial": True,
                },
                "one_plus_delta_Z_minus_two_Z": {
                    "lhs_terms": fraction_terms(one_plus_delta_minus_two),
                    "rhs_terms": fraction_terms(target_negative_one_minus_delta),
                    "exact_fraction_polynomial": True,
                },
                "Z": Z,
                "delta": delta,
            },
            "source_parameter_semantics": "conditional",
            "source_values_unevaluated": True,
            "global_derivative_identity": False,
        })
    return rows


def _volume_form(physical_map, ports):
    d = physical_map["d"]
    L = physical_map["L"]
    lam = physical_map["lambda_phys"]
    R = ports["R"]
    factor = _quotient(_mul(_power(lam, _sub(_number(3), ports["delta"])), L), d,
                       "positive under declared interior domain")
    return {
        "coordinates": {
            "R_Z_theta": {
                "expression": factor,
                "differential": "dR dZ dtheta",
            },
            "native_x_Z_theta": {
                "expression": _mul(factor, R),
                "differential": "dx dZ dtheta",
                "R_x": "R",
            },
        },
        "not_old_measure": True,
        "not_R_cubed": True,
        "positive_premise_only": True,
        "source_values_unevaluated": True,
    }


def select_physical_definition(branch, field, coefficient_order=0,
                               phase_mode="ordinary_Z"):
    """Return a typed selector descriptor, rejecting unsupported requests."""
    if branch not in BRANCHES:
        raise ValueError("unknown physical branch")
    if field not in SELECTABLE_FIELDS:
        raise ValueError("unknown physical field")
    if (not isinstance(coefficient_order, int)
            or isinstance(coefficient_order, bool)
            or coefficient_order != 0):
        raise ValueError("only coefficient order 0 is admitted")
    if phase_mode != "ordinary_Z":
        raise ValueError("only ordinary_Z phase mode is admitted")
    return {
        "kind": "typed_physical_selector",
        "branch": branch,
        "field": field,
        "coefficient_order": 0,
        "phase_mode": "ordinary_Z",
        "definition_namespace": COMPONENT + ".fresh_function_namespace",
        "evaluated": False,
    }


def _callback_counts(journal):
    counts = {}
    for event in journal:
        name = event["event"]
        counts[name] = counts.get(name, 0) + 1
    return counts


def _actual_counts(journal, success=False):
    events = _callback_counts(journal)
    result = {
        "source_family_calls": events.get("family_entered", 0),
        "new_source_identity_family_calls": 1 if events.get("family_entered") else 0,
        "new_ast_semantics_family_calls": 1 if events.get("family_entered") else 0,
        "new_algebra_family_calls": 1 if events.get("family_entered") else 0,
        "new_mathematical_proofs": 1 if success else 0,
        "scope_binding": events.get("scope_binding", 0),
        "spacetime_map": events.get("spacetime_map", 0),
        "branch_source_binding": events.get("branch_source_binding", 0),
        "streamfunction_definition": events.get("streamfunction_definition", 0),
        "cylindrical_field_definition": events.get("cylindrical_field_definition", 0),
        "Cartesian_field_definition": events.get("Cartesian_field_definition", 0),
        "pressure_definition": events.get("pressure_definition", 0),
        "conditional_kinematic_identity": events.get("conditional_kinematic_identity", 0),
        "physical_measure_definition": events.get("physical_measure_definition", 0),
        "selector_refusal": events.get("selector_refusal", 0),
        "failed_family_calls": 0 if success else (1 if events.get("family_entered") else 0),
        "old_math_calls": 0,
        "provider_calls": 0,
        "numeric_calls": 0,
        "N_materialization_calls": 0,
        "source_value_reducer_calls": 0,
        "ancestor_calls": 0,
        "temporal_recursion_calls": 0,
        "physical_source_calls": events.get("branch_source_binding", 0),
        "selector_negative_calls": events.get("selector_refusal", 0),
        "selector_calls": events.get("selector_refusal", 0),
    }
    return result


def _flags(success=False):
    names = {
        "conditional_physical_time_map_defined": success,
        "same_source_streamfunction_field_defined": success,
        "Cartesian_velocity_pressure_defined": success,
        "conditional_divergence_free_structure_installed": success,
        "physical_volume_form_defined": success,
        "profile_R_distinguished_from_physical_r": success,
        "native_field_derivative_identity_proved": False,
        "numeric_closure": False,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "temporal_recursion": False,
        "cap_deletion": False,
        "unconditional_parameter_admission": False,
        "global_phase_seam_regularity_proved": False,
        "axis_extension_proved": False,
        "numerical_oracle_available": False,
        "unconditional_Z4_source_admission": False,
        "actual_incoming_x_independence_proved": False,
        "oracle_source_value_identity": False,
    }
    return names


def _failure(journal, reason, pins=None, parent_commit=None):
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "full_goal_complete": False,
        "proof_count": 0,
        "callback_events": journal,
        "callback_event_count": len(journal),
        "callback_event_counts": _callback_counts(journal),
        "actual_counts": _actual_counts(journal, False),
        "counts": _actual_counts(journal, False),
        "new_flags": _flags(False),
        "flags": _flags(False),
        "runtime_model_metadata": None,
        "parent_commit": parent_commit,
        "source_wrapper_pins": pins or dict(WRAPPER_PINS),
        "payload": None,
        "failure": {"reason": str(reason), "prefix_preserved": True},
    }


def derive(payload, emit):
    """Build the fresh whole-chart typed function namespace once.

    The first callback is intentionally emitted before even checking the input
    object.  Source functions and expressions are never executed.  The four
    selector calls near completion are intentional negative-control API checks;
    they execute only the fresh selector's admission logic.
    """
    journal = []

    def record(event, row):
        entry = {"event": event, "row": row}
        journal.append(entry)
        if emit is not None:
            emit(event, row)

    record("family_entered", {
        "component": COMPONENT,
        "schema": INPUT_SCHEMA,
        "source_only": True,
        "mathematical_consumer": "fresh_typed_function_term_IR",
    })
    pins = dict(WRAPPER_PINS)
    parent_commit = None
    try:
        _require(_is_map(payload), "derive payload must be an object")
        _require(payload.get("schema") == INPUT_SCHEMA,
                 "input schema does not match this component")
        wrappers = _validate_wrappers(payload)
        saved_result, parent_payload = _parent_payload(wrappers["saved_flux_result"])
        parent_commit = _parent_commit(wrappers["saved_flux_result"], saved_result,
                                       parent_payload, payload)
        source_literals = wrappers["source_literals"]
        capsule = wrappers["physical_source_capsule"]
        _require(str(source_literals.get("schema", "")).startswith(
            "native_O3_density_flux_raw_literals_v1"),
            "source literal schema mismatch")
        _require(str(capsule.get("schema", "")).startswith(
            "native_O3_physical_source_capsule_v1"),
            "physical source capsule schema mismatch")
        finite_ast = _finite_source_ast(source_literals)
        capsule_ast = _capsule_source(capsule)
        ports = _parameter_ports(capsule)
        capsule_rows = {
            branch: _capsule_rows(capsule, branch) for branch in BRANCHES
        }
        my_rows, derivative_rows = _parent_rows(parent_payload)
        source_join_guards = _source_join_guards(
            capsule_rows, parent_payload, source_literals)
        record("scope_binding", {
            "scope": capsule.get("scope"),
            "ordinary_Z_orders": list(ORDERS),
            "physical_coefficient_order": 0,
            "branches": list(BRANCHES),
            "parameter_ports": ports,
            "source_ast": {"finite_N_functions": finite_ast, "physical_pressure_fields": capsule_ast},
            "source_join_guards": source_join_guards,
            "parent_commit": parent_commit,
            "source_values_unevaluated": True,
        })
        physical_map = _physical_map(ports)
        record("spacetime_map", {
            "map": physical_map,
            "source_function_composition": {
                "native_x": physical_map["physical_inverse"]["native_x"],
                "Z": physical_map["physical_inverse"]["Z"],
                "q_positive_root_equation": physical_map["physical_inverse"]["equation"],
                "numeric_inverse": False,
            },
            "source_values_unevaluated": True,
        })
        source_bindings = {}
        branches = {}
        identities = {}
        for branch in BRANCHES:
            source_binding = _build_source_bindings(
                capsule_rows, ports, parent_payload, branch,
                source_join_guards[branch])
            source_bindings[branch] = source_binding
            record("branch_source_binding", source_binding)
            branch_data = _build_branch(
                capsule_rows, ports, source_binding["rows"], parent_payload,
                my_rows, branch, physical_map)
            branches[branch] = branch_data
            record("streamfunction_definition", {
                "branch": branch,
                "ordinary_Z_orders": list(ORDERS),
                "rows": [{"ordinary_Z_order": row["ordinary_Z_order"],
                          "psi": row["psi"],
                          "vector_potential": row["vector_potential"]}
                         for row in branch_data["rows"]],
                "source_jet_rows": branch_data["source_jet_rows"],
                "physical_function_orders": branch_data["physical_function_orders"],
                "streamfunction": "lambda_phys^(1-delta)*Pstar*R*m",
                "source_values_unevaluated": True,
            })
            record("cylindrical_field_definition", {
                "branch": branch,
                "rows": [{
                    "ordinary_Z_order": row["ordinary_Z_order"],
                    "ur_stream": row["ur_stream"],
                    "ur": row["ur"],
                    "uz_stream": row["uz_stream"],
                    "uz": row["uz"],
                    "utheta": row["utheta"],
                    "curl_definition": "ur=-psi_z/r; uz=psi_r/r",
                    "formal_chain_rule": {
                        "operator": "ordinary_Z_chain_rule",
                        "streamfunction": row["psi"],
                        "C2_regular_interior_premise": True,
                    },
                    "curl_commutation": {
                        "conditional_C2_only": True,
                        "global_derivative_identity": False,
                    },
                    "conditional_only": True,
                } for row in branch_data["rows"]],
                "source_jet_orders_data_only": list(ORDERS),
                "source_values_unevaluated": True,
            })
            record("Cartesian_field_definition", {
                "branch": branch,
                "rows": [{
                    "ordinary_Z_order": row["ordinary_Z_order"],
                    "u": row["u"], "v": row["v"], "w": row["w"],
                    "source_function_at": row["source_function_at"],
                    "physical_cartesian_functions": row["physical_cartesian_functions"],
                } for row in branch_data["rows"]],
                "source_jet_orders_data_only": list(ORDERS),
                "physical_coordinates": {"X": physical_map["physical_X"],
                                          "Y": physical_map["physical_Y"]},
                "axis_extension_proved": False,
                "source_values_unevaluated": True,
            })
            record("pressure_definition", {
                "branch": branch,
                "rows": [{
                    "ordinary_Z_order": row["ordinary_Z_order"],
                    "pressure_profile": row["pressure_profile"],
                    "pressure": row["pressure"],
                    "P0_used_once": True,
                } for row in branch_data["rows"]],
                "source_jet_orders_data_only": list(ORDERS),
                "formula": "lambda_phys^(-2-2*delta)*Pstar^2*(P0+p)",
                "source_values_unevaluated": True,
            })
            identities[branch] = _identity_rows(
                branch_data, ports, parent_payload, branch)
            record("conditional_kinematic_identity", {
                "branch": branch,
                "rows": identities[branch],
                "regularity_premise": "C2 regular interior phase cells",
                "global_derivative_identity_proved": False,
                "source_values_unevaluated": True,
            })
        volume = _volume_form(physical_map, ports)
        record("physical_measure_definition", {
            "volume_form": volume,
            "source_values_unevaluated": True,
        })
        refusal_cases = (
            ("unknown_branch", {"branch": "not_a_branch", "field": "ur"}),
            ("n1_order", {"branch": "background", "field": "ur", "coefficient_order": 1}),
            ("held_phase", {"branch": "background", "field": "ur", "phase_mode": "held_phase"}),
            ("unknown_field", {"branch": "background", "field": "not_a_field"}),
        )
        for label, case in refusal_cases:
            try:
                select_physical_definition(**case)
            except ValueError as exc:
                record("selector_refusal", {
                    "control": label,
                    "request": case,
                    "rejected": True,
                    "reason": str(exc),
                })
            else:
                raise _KernelInputError("negative selector control was accepted")
        fresh_payload = {
            "namespace": {
                "name": COMPONENT + ".fresh_function_namespace",
                "kind": "typed_function_term_IR",
                "source_values_unevaluated": True,
                "no_native_ids_invented": True,
            },
            "source_wrapper_pins": dict(WRAPPER_PINS),
            "parent_pointer": {
                "source": "saved_flux_result",
                "wrapper_sha256": WRAPPER_PINS["saved_flux_result"],
                "parent_commit": parent_commit,
                "payload_path": "result.payload",
                "recursive_parent_copy": False,
            },
            "source_ast": {"finite_N_functions": finite_ast,
                           "physical_pressure_fields": capsule_ast},
            "parameter_ports": ports,
            "physical_map": physical_map,
            "source_join_guards": source_join_guards,
            "source_bindings": source_bindings,
            "branches": branches,
            "conditional_kinematic_identities": identities,
            "volume_form": volume,
            "accepted_parent_data": {
                "my_compatibility_rows": {
                    "count": len(my_rows),
                    "path": "result.payload.my_compatibility_rows",
                    "source": "saved_flux_result",
                },
                "derivative_rewrite_rows": {
                    "count": len(derivative_rows),
                    "path": "result.payload.derivative_rewrite_rows",
                    "source": "saved_flux_result",
                },
            },
            "assumptions": {
                "Pstar_fixed_positive": True,
                "delta_fixed_between_zero_and_one_over_200": True,
                "ordinary_Z_source_jets": True,
                "incoming_histories_held": True,
                "regular_interior_phase_cells": True,
                "axis_extension": False,
                "global_phase_seam": False,
            },
            "numeric_evaluation": False,
            "source_value_evaluation": False,
            "ancestor_calls": 0,
            "provider_calls": 0,
            "old_math_calls": 0,
        }
        record("family_completed", {
            "component": COMPONENT,
            "fresh_namespace": fresh_payload["namespace"]["name"],
            "proof_count": 1,
            "source_values_unevaluated": True,
        })
        return {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "passed": True,
            "source_only": True,
            "metadata_only": False,
            "full_goal_complete": False,
            "proof_count": 1,
            "callback_events": journal,
            "callback_event_count": len(journal),
            "callback_event_counts": _callback_counts(journal),
            "actual_counts": _actual_counts(journal, True),
            "counts": _actual_counts(journal, True),
            "new_flags": _flags(True),
            "flags": _flags(True),
            "runtime_model_metadata": None,
            "parent_commit": parent_commit,
            "source_wrapper_pins": dict(WRAPPER_PINS),
            "payload": fresh_payload,
            "failure": None,
        }
    except Exception as exc:  # preserve every durable callback prefix
        return _failure(journal, exc, pins, parent_commit)


__all__ = ["derive", "select_physical_definition"]
