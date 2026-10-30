"""Exact conditional Rc recipe comparison extracted from pinned source ASTs.

No archived graph is read and no historical module is imported or called.
Source-level role continuity and common parameter ownership remain explicit
premises; the symbolic identity is not a physical or endpoint certificate.
"""
from __future__ import annotations

import ast
import sympy as s


class UnsupportedSource(ValueError):
    pass


def function(tree, name, class_name=None):
    body = tree.body
    if class_name:
        matches = [n for n in body if isinstance(n, ast.ClassDef) and n.name == class_name]
        if len(matches) != 1:
            raise UnsupportedSource("missing/ambiguous class " + class_name)
        body = matches[0].body
    found = [n for n in body if isinstance(n, ast.FunctionDef) and n.name == name]
    if len(found) != 1:
        raise UnsupportedSource("missing/ambiguous function " + name)
    return found[0]


def assignments(fn, target):
    return [node for node in ast.walk(fn) if isinstance(node, ast.Assign)
            and any(ast.unparse(t) == target for t in node.targets)]


def unique_assignment(fn, target):
    found = assignments(fn, target)
    if len(found) != 1:
        raise UnsupportedSource("missing/ambiguous assignment " + target)
    return found[0].value


def named_dict_value(node, key):
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "dict":
        values = [kw.value for kw in node.keywords if kw.arg == key]
    elif isinstance(node, ast.Dict):
        values = [value for k, value in zip(node.keys, node.values)
                  if isinstance(k, ast.Constant) and k.value == key]
    else:
        raise UnsupportedSource("literal source dictionary required")
    if len(values) != 1:
        raise UnsupportedSource("missing/ambiguous literal entry " + key)
    return values[0]


def evaluate(node, env):
    """Closed arithmetic-only interpreter; unknown calls and symbols fail."""
    if isinstance(node, ast.Constant) and not isinstance(node.value, bool):
        if isinstance(node.value, (int, str)):
            return s.Rational(node.value)
        raise UnsupportedSource("inexact literal")
    key = ast.unparse(node)
    if key in env:
        return env[key]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -evaluate(node.operand, env)
    if isinstance(node, ast.BinOp):
        a, b = evaluate(node.left, env), evaluate(node.right, env)
        if isinstance(node.op, ast.Add): return a + b
        if isinstance(node.op, ast.Sub): return a - b
        if isinstance(node.op, ast.Mult): return a * b
        if isinstance(node.op, ast.Div): return a / b
        raise UnsupportedSource("unsupported binary operator")
    if isinstance(node, ast.Call) and not node.keywords:
        op = ast.unparse(node.func)
        if op == "g.unary" and len(node.args) == 2:
            if not isinstance(node.args[0], ast.Constant): raise UnsupportedSource("unary name must be literal")
            return evaluate(ast.Call(func=ast.Name(id="sy." + str(node.args[0].value)), args=[node.args[1]], keywords=[]), env)
        args = [evaluate(v, env) for v in node.args]
        if op in ("c", "g.constant", "sy.Integer") and len(args) == 1: return args[0]
        if op == "sy.Rational" and len(args) == 2: return args[0] / args[1]
        if op in ("sy.log", "g.log") and len(args) == 1: return s.log(args[0])
        if op in ("sy.exp", "g.exp") and len(args) == 1: return s.exp(args[0])
        if op == "g.add": return s.Add(*args)
        if op == "g.mul": return s.Mul(*args)
        if op == "g.sub" and len(args) == 2: return args[0] - args[1]
        if op == "g.neg" and len(args) == 1: return -args[0]
    raise UnsupportedSource("unsupported expression: " + ast.unparse(node))


def identical_ast(actual, expected):
    return ast.dump(actual, include_attributes=False) == ast.dump(ast.parse(expected, mode="eval").body, include_attributes=False)


def prove(geometry_text, frame_text, prescribed_text, outer_text):
    geometry, frame, prescribed, outer = map(ast.parse, (geometry_text, frame_text, prescribed_text, outer_text))
    C, P, T, B, S, sc, x, p = s.symbols("logC logP T hbB hbS sc x p", real=True)
    env = {"C": C, "P": P, "T": T, "B": B, "S": S, "sc": sc, "x": x}
    geom_fn = function(geometry, "exact_geometry")
    for name in ("Y", "L", "Rm"):
        env[name] = evaluate(unique_assignment(geom_fn, name), env)
    old_offset = evaluate(named_dict_value(unique_assignment(geom_fn, "offsets"), "O3_power"), env)
    updates = [n for n in ast.walk(geom_fn) if isinstance(n, ast.Call) and ast.unparse(n.func) == "domains.update"]
    o3_domain = [kw.value for n in updates for kw in n.keywords if kw.arg == "O3_power"]
    if len(o3_domain) != 1 or ast.unparse(o3_domain[0]) != "(sy.Integer(0), sy.Integer(2))":
        raise UnsupportedSource("old O3 domain is not [0,2]")

    frame_fn = function(frame, "build")
    old_env = {"params['logP']": P, "params['hbB']": B, "params['sc']": sc}
    logminus = evaluate(unique_assignment(frame_fn, "logminus"), old_env)
    if not identical_ast(unique_assignment(frame_fn, "rc"), "g.add(phase_rc, logminus)"):
        raise UnsupportedSource("old Rc composition changed")
    old_logrc = old_offset.subs(x, 2) + logminus

    owner = function(prescribed, "__init__", "CurrentOriginalPrescribedR100RmC6Functions")
    initial = {"self.core.parameters['logC']": C, "self.core.parameters['logP']": P, "p": p}
    logref = evaluate(unique_assignment(owner, "logref"), initial)
    offsets = [a for a in assignments(owner, "offset") if identical_ast(a.value, "g.add(c(-7), p)")]
    radii = [a for a in assignments(owner, "R") if identical_ast(a.value, "g.mul(c(110), g.unary('exp', g.add(logref, offset)))")]
    if len(offsets) != 1 or len(radii) != 2:
        raise UnsupportedSource("postrestore source shape changed")
    post_offset = evaluate(offsets[0].value, initial)
    rm_radius = evaluate(radii[-1].value, {"logref": logref, "offset": post_offset}).subs(p, 1)
    rm_assignments = [n for n in ast.walk(owner) if isinstance(n, ast.Assign)
                      and any(ast.unparse(t) == "self.endpoints['Rm']" for t in n.targets)]
    if len(rm_assignments) != 1 or not identical_ast(rm_assignments[0].value, "self.end('post_restore', g.one)"):
        raise UnsupportedSource("Rm endpoint role changed")

    outer_geom = function(outer, "geometry")
    current_Y = evaluate(unique_assignment(outer_geom, "Y"), {})
    offset_subscript = unique_assignment(outer_geom, "offset")
    if not isinstance(offset_subscript, ast.Subscript): raise UnsupportedSource("outer geometry dispatch changed")
    outer_offset = evaluate(named_dict_value(offset_subscript.value, "O3_power"), {"Y": current_Y, "x": x})
    domains = unique_assignment(outer_geom, "domains")
    if ast.unparse(named_dict_value(domains, "O3_power")) != "(g.zero, c(2))":
        raise UnsupportedSource("current O3 domain is not [0,2]")
    outer_owner = function(outer, "__init__", "CurrentOriginalTrueOuterC6C5Functions")
    if not identical_ast(unique_assignment(outer_owner, "self.Rm"), "self.leading.inlet['radius']"):
        raise UnsupportedSource("current Rm owner changed")
    outer_r = unique_assignment(outer_owner, "R")
    if not identical_ast(outer_r, "g.mul(self.Rm, g.unary('exp', offset))"):
        raise UnsupportedSource("current outer radius expression changed")
    current_rc = rm_radius * s.exp(outer_offset.subs(x, 2))
    log_current_rc = s.log(110) + logref + post_offset.subs(p, 1) + outer_offset.subs(x, 2)
    residual = s.simplify(s.expand_log(old_logrc - log_current_rc, force=True))
    if residual != 0:
        raise UnsupportedSource("old/current Rc source recipes differ")
    return {"status": "conditional_formula_identity_proved",
            "old_logRc": str(s.simplify(old_logrc)), "current_logRc": str(s.simplify(log_current_rc)),
            "difference": str(residual), "symbolic_domain": "real logC, logP, hbB, sc; positive rational logarithm arguments",
            "old_O3_endpoint": 2, "current_postrestore_endpoint": 1, "current_O3_endpoint": 2,
            "source_assumptions": ["old phase_rc is the exact O3 right offset", "old/current logC and logP have the same source ownership and units", "current leading inlet uses the prescribed Rm endpoint", "endpoint substitutions have their defining function semantics"],
            "micro_origin_terms_retained_and_cancel_in_recipe": True,
            "current_rc_positive_from_recipe": bool(current_rc.is_positive),
            "cross_owner_source_identity_proved": False, "physical_units_proved": False,
            "endpoint_identity_proved": False, "numeric_source_values_computed": False,
            "old_errors_removed": False, "global_closure_proved": False}
