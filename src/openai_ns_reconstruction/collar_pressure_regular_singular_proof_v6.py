"""Conditional source algebra for the original collar pressure seam.

The public API accepts caller-pinned source text and performs no file I/O.  It
uses the actual AST assignments in the two fixed producer sources, then uses
SymPy only for the conditional algebra on their common collar domain.  It does
not assert cross-owner source identity, a numerical value, equality of the
archived finite six rows, or error cancellation.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Any, Callable, Mapping

import sympy as sp


FORWARD_SOURCE_PATH = "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py"
CANONICAL_SOURCE_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py"
)
SOURCE_PATHS = (FORWARD_SOURCE_PATH, CANONICAL_SOURCE_PATH)


class SourceExtractionError(ValueError):
    """Raised when the pinned AST cannot be extracted without guessing."""


@dataclass(frozen=True)
class _Field:
    value: sp.Expr
    Z: sp.Expr
    ZZ: sp.Expr


@dataclass(frozen=True)
class _Namespace:
    name: str


@dataclass(frozen=True)
class _Method:
    owner: str
    name: str


def _num(value: Any) -> sp.Expr:
    if isinstance(value, sp.Expr):
        return value
    if isinstance(value, bool):
        return sp.Integer(int(value))
    if isinstance(value, int):
        return sp.Integer(value)
    if isinstance(value, float):
        return sp.Float(value)
    if isinstance(value, str):
        try:
            return sp.Rational(value)
        except (TypeError, ValueError):
            raise SourceExtractionError(f"non-numeric source string: {value!r}") from None
    raise SourceExtractionError(f"unsupported source constant: {value!r}")


def _expr(value: Any) -> sp.Expr:
    return value.value if isinstance(value, _Field) else _num(value)


def _parse(text: str, filename: str) -> ast.Module:
    if not isinstance(text, str) or not text.strip():
        raise SourceExtractionError(f"empty pinned source: {filename}")
    try:
        return ast.parse(text, filename=filename, mode="exec")
    except (SyntaxError, TypeError) as exc:
        raise SourceExtractionError(f"cannot parse {filename}: {exc}") from exc


def _function(tree: ast.Module, name: str, filename: str) -> ast.FunctionDef:
    found = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    if len(found) != 1:
        raise SourceExtractionError(f"{filename}: expected one top-level def {name}")
    return found[0]


def _assignment(function: ast.FunctionDef, name: str) -> ast.Assign | ast.AnnAssign:
    found = []
    for node in function.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(part, ast.Name) and part.id == name
            for target in node.targets for part in ast.walk(target)
        ):
            found.append(node)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == name:
            found.append(node)
    if len(found) != 1:
        raise SourceExtractionError(f"{function.name}: expected one assignment {name}")
    return found[0]


def _value(node: ast.Assign | ast.AnnAssign) -> ast.AST:
    if isinstance(node, ast.Assign):
        if len(node.targets) != 1:
            raise SourceExtractionError("chained assignment is not accepted")
        return node.value
    if node.value is None:
        raise SourceExtractionError("annotation without a value is not accepted")
    return node.value


def _slice(node: ast.Subscript) -> ast.AST:
    value = node.slice
    if isinstance(value, ast.Index):  # pragma: no cover - legacy parser
        return value.value
    return value


class _Eval:
    """Small, closed evaluator for the arithmetic DSL used by the sources."""

    def __init__(
        self,
        environment: Mapping[str, Any],
        z: sp.Symbol,
        q_operator: Callable[[sp.Expr, sp.Expr, sp.Expr], sp.Expr] | None = None,
    ) -> None:
        self.env = dict(environment)
        self.z = z
        self.q_operator = q_operator

    def run(self, node: ast.AST, env: Mapping[str, Any] | None = None) -> Any:
        e = self.env if env is None else env
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in e:
                raise SourceExtractionError(f"unbound source name: {node.id}")
            return e[node.id]
        if isinstance(node, ast.List):
            result = []
            for item in node.elts:
                if isinstance(item, ast.Starred):
                    result.extend(self.run(item.value, e))
                else:
                    result.append(self.run(item, e))
            return result
        if isinstance(node, ast.Tuple):
            return tuple(self.run(item, e) for item in node.elts)
        if isinstance(node, ast.Dict):
            if any(key is None for key in node.keys):
                raise SourceExtractionError("dictionary unpacking is not accepted")
            return {self.run(k, e): self.run(v, e) for k, v in zip(node.keys, node.values)}
        if isinstance(node, ast.Slice):
            return slice(*(self.run(v, e) if v is not None else None
                           for v in (node.lower, node.upper, node.step)))
        if isinstance(node, ast.Subscript):
            base, key = self.run(node.value, e), self.run(_slice(node), e)
            try:
                return base[key]
            except (KeyError, IndexError, TypeError) as exc:
                raise SourceExtractionError(f"unresolved subscript: {ast.unparse(node)}") from exc
        if isinstance(node, ast.Attribute):
            base = self.run(node.value, e)
            if isinstance(base, _Namespace):
                if node.attr == "one":
                    return sp.Integer(1)
                if node.attr == "zero":
                    return sp.Integer(0)
                return _Method(base.name, node.attr)
            if isinstance(base, _Field):
                if node.attr in {"value", "Z", "ZZ"}:
                    return getattr(base, node.attr)
                if node.attr == "dz":
                    return lambda: base.Z
            if isinstance(base, sp.Expr):
                if node.attr in {"value", "Z"}:
                    return base if node.attr == "value" else sp.diff(base, self.z)
                if node.attr == "ZZ":
                    return sp.diff(base, self.z, 2)
                if node.attr == "dz":
                    return lambda: sp.diff(base, self.z)
            raise SourceExtractionError(f"unsupported source attribute: {ast.unparse(node)}")
        if isinstance(node, ast.Lambda):
            names = [arg.arg for arg in node.args.args]

            def closure(*args: Any) -> Any:
                if len(args) != len(names):
                    raise SourceExtractionError("source lambda arity changed")
                child = dict(e)
                child.update(zip(names, args))
                return self.run(node.body, child)

            return closure
        if isinstance(node, ast.BinOp):
            left, right = _expr(self.run(node.left, e)), _expr(self.run(node.right, e))
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                return left**right
            raise SourceExtractionError("unsupported arithmetic operator")
        if isinstance(node, ast.UnaryOp):
            value = _expr(self.run(node.operand, e))
            if isinstance(node.op, ast.USub):
                return -value
            if isinstance(node.op, ast.UAdd):
                return value
            raise SourceExtractionError("unsupported unary operator")
        if isinstance(node, ast.BoolOp):
            if isinstance(node.op, ast.Or):
                return next(value for item in node.values if (value := self.run(item, e)))
            if isinstance(node.op, ast.And):
                value = None
                for item in node.values:
                    value = self.run(item, e)
                    if not value:
                        return value
                return value
            raise SourceExtractionError("unsupported boolean operator")
        if isinstance(node, ast.Call):
            fn = self.run(node.func, e)
            args = [self.run(arg, e) for arg in node.args]
            kwargs = {
                keyword.arg: self.run(keyword.value, e)
                for keyword in node.keywords
                if keyword.arg is not None
            }
            if isinstance(fn, _Method):
                return self._method(fn.owner, fn.name, args, kwargs)
            if not callable(fn):
                raise SourceExtractionError(f"non-callable source expression: {ast.unparse(node)}")
            try:
                return fn(*args, **kwargs)
            except SourceExtractionError:
                raise
            except Exception as exc:
                raise SourceExtractionError(f"source call failed: {ast.unparse(node)}: {exc}") from exc
        if isinstance(node, ast.ListComp):
            return self._comprehension(node.elt, node.generators, e)
        if isinstance(node, ast.GeneratorExp):
            return tuple(self._comprehension(node.elt, node.generators, e))
        raise SourceExtractionError(f"unsupported source AST node: {type(node).__name__}")

    def _comprehension(self, elt: ast.AST, generators: list[ast.comprehension], env: Mapping[str, Any]) -> list[Any]:
        if len(generators) != 1 or generators[0].is_async:
            raise SourceExtractionError("only one synchronous comprehension is accepted")
        generator = generators[0]
        result = []
        for item in self.run(generator.iter, env):
            child = dict(env)
            target = generator.target
            if isinstance(target, ast.Name):
                child[target.id] = item
            elif isinstance(target, ast.Tuple) and isinstance(item, (tuple, list)):
                if len(target.elts) != len(item) or not all(isinstance(x, ast.Name) for x in target.elts):
                    raise SourceExtractionError("comprehension destructuring changed")
                child.update({x.id: value for x, value in zip(target.elts, item)})
            else:
                raise SourceExtractionError("comprehension target changed")
            if all(self.run(test, child) for test in generator.ifs):
                result.append(self.run(elt, child))
        return result

    def _method(self, owner: str, name: str, args: list[Any], kwargs: Mapping[str, Any]) -> Any:
        if owner in {"a", "g"}:
            if owner == "g" and name == "unary":
                if len(args) != 2 or not isinstance(args[0], str):
                    raise SourceExtractionError("g.unary signature changed")
                value = _expr(args[1])
                return {"exp": sp.exp, "log": sp.log, "sqrt": sp.sqrt}.get(args[0], lambda _: None)(value)
            if owner == "g" and name == "constant":
                if len(args) != 1:
                    raise SourceExtractionError("g.constant signature changed")
                return _num(args[0])
            if owner == "g" and name == "symbol":
                if len(args) != 1 or not isinstance(args[0], str):
                    raise SourceExtractionError("g.symbol signature changed")
                return self.env.get(args[0], sp.Symbol(args[0]))
            if owner == "g" and name == "node":
                if kwargs.get("name") == "original_flat_sigma_prime":
                    return self.env["sigma_prime"]
                raise SourceExtractionError("unexpected g.node source")
            if name in {"div", "quotient"} and len(args) == 3:
                if not isinstance(args[2], str):
                    raise SourceExtractionError("quotient reason must be a literal string")
                args = args[:2]
            values = [_expr(value) for value in args]
            if name == "add":
                return sum(values, sp.Integer(0))
            if name == "mul":
                return sp.prod(values)
            if name in {"sub", "div", "quotient"} and len(values) >= 2:
                return values[0] - values[1] if name == "sub" else values[0] / values[1]
            if name == "scale" and len(values) == 2:
                return values[0] * values[1]
            if name == "fixed" and len(values) == 1:
                return values[0]
            if name == "neg" and len(values) == 1:
                return -values[0]
            raise SourceExtractionError(f"unsupported {owner}.{name} source method")
        if owner == "source" and name == "C1Function":
            if not args:
                raise SourceExtractionError("empty C1Function")
            return _expr(args[0])
        if owner == "parent" and name == "C1Algebra":
            return _Namespace("a")
        if owner == "c3" and name == "C3Algebra":
            return _Namespace("a")
        if owner == "own" and name == "radial_Q":
            if self.q_operator is None or len(args) != 4:
                raise SourceExtractionError("canonical radial_Q operator unavailable")
            return self.q_operator(_expr(args[1]), _expr(args[2]), _expr(args[3]))
        raise SourceExtractionError(f"unsupported {owner}.{name} source method")


def _set_target(node: ast.Assign | ast.AnnAssign, value: Any, env: dict[str, Any]) -> None:
    target = node.target if isinstance(node, ast.AnnAssign) else node.targets[0]
    if isinstance(target, ast.Name):
        env[target.id] = value
        return
    if isinstance(target, (ast.Tuple, ast.List)) and isinstance(value, (tuple, list)):
        if len(target.elts) != len(value) or not all(isinstance(x, ast.Name) for x in target.elts):
            raise SourceExtractionError("assignment destructuring changed")
        env.update({x.id: item for x, item in zip(target.elts, value)})
        return
    raise SourceExtractionError("unsupported assignment target")


def _until(function: ast.FunctionDef, evaluator: _Eval, env: dict[str, Any], stop: str) -> dict[str, Any]:
    for node in function.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            _set_target(node, evaluator.run(_value(node), env), env)
            if stop in env:
                return env
            continue
        if isinstance(node, ast.Return):
            raise SourceExtractionError(f"{function.name}: return before {stop}")
        raise SourceExtractionError(f"{function.name}: unsupported statement before {stop}")
    raise SourceExtractionError(f"{function.name}: missing {stop}")


def _hit(function: ast.FunctionDef, name: str, expression: Any) -> dict[str, Any]:
    node = _assignment(function, name)
    return {
        "function": function.name,
        "assignment": name,
        "line": node.lineno,
        "ast": ast.dump(_value(node), include_attributes=False),
        "formula": ("callable source expression" if callable(expression) else
                    sp.sstr(_expr(expression)) if not isinstance(expression, (list, tuple)) else [sp.sstr(_expr(x)) for x in expression]),
    }


def _normalize(value: Any, raw_rho: sp.Expr, rho: sp.Symbol) -> Any:
    if isinstance(value, list):
        return [_normalize(item, raw_rho, rho) for item in value]
    if isinstance(value, tuple):
        return tuple(_normalize(item, raw_rho, rho) for item in value)
    if isinstance(value, sp.Expr):
        return sp.simplify(value.subs(raw_rho / 4, rho / 4).xreplace({raw_rho: rho}))
    return value


def _forward(text: str) -> tuple[dict[str, Any], dict[str, Any], sp.Expr]:
    tree = _parse(text, FORWARD_SOURCE_PATH)
    leading = _function(tree, "leading_regular_frame", FORWARD_SOURCE_PATH)
    pressure = _function(tree, "normalized_pressure_source", FORWARD_SOURCE_PATH)
    z, rho = sp.symbols("Z rho", commutative=True)
    eps, hb, de, P2 = sp.symbols("epsilon hb delta Pstar2", nonzero=True)
    sigma, sigma_prime = sp.symbols("sigma sigma_prime")
    phi0, phi1, phi2 = sp.symbols("Phi_core Phi_core_rho Phi_core_rhorho")
    psi0, psi1, psi2 = sp.symbols("Psi_actual Psi_core_rho Psi_core_rhorho")
    mean0, ratio, Ubase, axis_j = sp.symbols("mean_actual ratio U_base axis_j")
    env = {
        "z": z, "rho": rho, "eps": eps, "de": de, "P2": P2, "hb": hb,
        "sigma": sigma, "sigma_prime": sigma_prime, "axis_j": axis_j,
        "g": _Namespace("g"), "a": _Namespace("a"), "c": _num,
        "source": _Namespace("source"), "parent": _Namespace("parent"),
        "own": _Namespace("own"), "controls": _Namespace("controls"),
        "c3": _Namespace("c3"), "zip": zip,
        "Laurent": lambda value: value[0] if isinstance(value, dict) and 0 in value else (_ for _ in ()).throw(SourceExtractionError("Laurent row changed")),
        "values": {"Phi": ratio * phi0, "Psi": psi0, "Mpsi": mean0},
        "comparison": {"Phi": [phi0, phi1, phi2], "Psi": [psi0, psi1, psi2]},
    }
    ev = _Eval(env, z)
    lead = _until(leading, ev, env, "Q")
    required = ("chi", "A", "Ar", "ratio", "phi", "pr", "psi", "mr", "mean", "U", "M", "d", "L", "Q")
    if any(name not in lead for name in required):
        raise SourceExtractionError("leading_regular_frame formula set is incomplete")
    q_node = _value(_assignment(leading, "Q"))
    if not isinstance(q_node, ast.ListComp) or len(q_node.generators) != 1:
        raise SourceExtractionError("forward Q is not the expected comprehension")
    generator = q_node.generators[0]
    if not (isinstance(generator.iter, ast.Call) and isinstance(generator.iter.func, ast.Name) and generator.iter.func.id == "zip" and isinstance(generator.target, ast.Tuple) and len(generator.target.elts) == 2):
        raise SourceExtractionError("forward Q generator changed")
    u_name, m_name = generator.target.elts
    if not isinstance(u_name, ast.Name) or not isinstance(m_name, ast.Name):
        raise SourceExtractionError("forward Q generator names changed")
    u_sym, m_sym, mz_sym = sp.symbols("Q_template_u Q_template_m Q_template_mz")
    qenv = dict(env)
    qenv.update({u_name.id: _Field(u_sym, u_sym, u_sym), m_name.id: _Field(m_sym, mz_sym, mz_sym)})
    q_template = _expr(ev.run(q_node.elt, qenv))

    def q_operator(de_expr: sp.Expr, u_expr: sp.Expr, m_expr: sp.Expr) -> sp.Expr:
        return q_template.subs({de: de_expr, u_sym: u_expr, m_sym: m_expr, mz_sym: sp.diff(m_expr, z)}, simultaneous=True)

    p_env = dict(env)
    p_env["frame"] = {"z": z, "rho": rho, "de": de, "e": eps, "P2": P2, "Q": lead["Q"], "U": lead["U"]}
    p_ev = _Eval(p_env, z)
    for statement in pressure.body:
        if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Constant):
            continue
        if isinstance(statement, (ast.Assign, ast.AnnAssign)):
            _set_target(statement, p_ev.run(_value(statement), p_env), p_env)
        elif isinstance(statement, ast.AugAssign) and isinstance(statement.target, ast.Name):
            expr = ast.BinOp(left=ast.Name(id=statement.target.id, ctx=ast.Load()),
                             op=statement.op, right=statement.value)
            p_env[statement.target.id] = p_ev.run(expr, p_env)
        elif isinstance(statement, ast.Return):
            break
        else:
            raise SourceExtractionError("unsupported normalized pressure statement")
    returns = [node for node in pressure.body if isinstance(node, ast.Return)]
    if len(returns) != 1:
        raise SourceExtractionError("normalized_pressure_source return changed")
    normalized = _expr(p_ev.run(returns[0].value, p_env))
    hits = {f"forward.{name}": _hit(leading, name, lead[name]) for name in required}
    hits["forward.normalized_pressure_source.return"] = {
        "function": pressure.name, "assignment": "return", "line": returns[0].lineno,
        "ast": ast.dump(returns[0].value, include_attributes=False), "formula": sp.sstr(normalized),
    }
    hits["forward.Q_template"] = {"function": leading.name, "assignment": "Q_template", "line": _assignment(leading, "Q").lineno, "ast": ast.dump(q_node.elt, include_attributes=False), "formula": sp.sstr(q_template)}
    payload = {
        "z": z, "rho": rho, "eps": eps, "hb": hb, "de": de, "P2": P2,
        "sigma": sigma, "sigma_prime": sigma_prime, "phi1": phi1, "psi1": psi1,
        "ratio": ratio, "leading": lead, "normalized": normalized,
        "normalized_physical": sp.simplify(P2 * normalized), "q_operator": q_operator,
        "Ubase": 4 * z + axis_j, "mean0": mean0,
    }
    return payload, hits, q_template


def _canonical(text: str, f: Mapping[str, Any], q_template: sp.Expr) -> tuple[dict[str, Any], dict[str, Any]]:
    tree = _parse(text, CANONICAL_SOURCE_PATH)
    collar = _function(tree, "collar_fields", CANONICAL_SOURCE_PATH)
    pressure = _function(tree, "direct_pressure_functions", CANONICAL_SOURCE_PATH)
    z, rho = f["z"], f["rho"]
    eps, hb, de, P2 = f["eps"], f["hb"], f["de"], f["P2"]
    sigma, sigma_prime, ratio = f["sigma"], f["sigma_prime"], f["ratio"]
    phi0, phi1, phi2 = sp.symbols("Phi_core Phi_core_rho Phi_core_rhorho")
    psi1, psi2 = f["psi1"], sp.symbols("Psi_core_rhorho")
    psi_actual, mean_actual = sp.symbols("Psi_actual mean_actual")
    phi4, psi4, mean4 = sp.symbols("Phi_core_at_4 Psi_core_at_4 mean_core_at_4")
    s = sp.Symbol("original_lower_collar_s")
    raw_rho = 4 * sp.exp(hb * s)
    core_rows = {"Phi": [phi0, phi1, phi2], "Psi": [psi4, psi1, psi2], "mean_Psi": [mean4, sp.Symbol("mean_core_rho"), sp.Symbol("mean_core_rhorho")]}
    parameters = {"hb": hb, "epsilon": eps, "delta": de, "sc_quarter": sp.Symbol("sc_quarter"), "sc_half": sp.Symbol("sc_half"), "logLambda": sp.Symbol("logLambda")}
    core_functions = {"U0": f["Ubase"], "F0": sp.Symbol("F0"), "P0": sp.Symbol("P0")}
    integral_symbols = {"logratio": sp.Symbol("logratio")}

    def at(_g: Any, expression: Any, _variable: Any, _value: Any) -> Any:
        return expression

    def integral(*_args: Any, **kwargs: Any) -> sp.Expr:
        if kwargs.get("original_ratio_ODE"):
            return integral_symbols["logratio"]
        if kwargs.get("actual_Psi_radial_transport"):
            return psi_actual - psi4
        if kwargs.get("actual_own_mean_FTC"):
            return rho * mean_actual - 4 * mean4
        return sp.Symbol("unresolved_history_integral")

    def exp3(_g: Any, _q: Any) -> sp.Expr:
        return ratio

    env = {
        "rho_variable": rho, "core_rows": core_rows, "core_functions": core_functions, "parameters": parameters,
        "g": _Namespace("g"), "a": _Namespace("a"), "c": _num, "source": _Namespace("source"),
        "parent": _Namespace("parent"), "own": _Namespace("own"), "controls": _Namespace("controls"), "c3": _Namespace("c3"), "zip": zip,
        "sigma": lambda _x: sigma, "sigma_prime": lambda _x: sigma_prime, "at": at, "integral": integral, "exp3": exp3,
    }
    qop = f["q_operator"]
    ev = _Eval(env, z, qop)
    state = _until(collar, ev, env, "Q")
    required = ("rho", "chi", "rs", "R", "jacobian", "ratio_body", "ratio", "Xr", "phi", "psi_body", "psi0", "psi", "mean0", "mean1", "mean2", "mean", "U", "m", "Q")
    if any(name not in state for name in required):
        raise SourceExtractionError("collar_fields formula set is incomplete")
    q_node = _value(_assignment(collar, "Q"))
    if not (isinstance(q_node, ast.ListComp) and len(q_node.generators) == 1 and isinstance(q_node.elt, ast.Call) and isinstance(q_node.elt.func, ast.Attribute) and isinstance(q_node.elt.func.value, ast.Name) and q_node.elt.func.value.id == "own" and q_node.elt.func.attr == "radial_Q" and len(q_node.elt.args) == 4):
        raise SourceExtractionError("canonical radial_Q delegation changed")
    state = {key: _normalize(value, raw_rho, rho) for key, value in state.items()}
    U0, Q0, Q1, Q2 = state["U"][0], state["Q"][0], state["Q"][1], state["Q"][2]
    fields = {
        "radius": eps * rho,
        "U": _Field(U0, sp.diff(U0, z), sp.diff(U0, z, 2)),
        "Q": _Field(Q0, sp.diff(Q0, z), sp.diff(Q0, z, 2)),
        "Q_y": _Field(rho * Q1, sp.diff(rho * Q1, z), sp.diff(rho * Q1, z, 2)),
        "Q_yy": _Field(rho**2 * Q2 + rho * Q1, sp.diff(rho**2 * Q2 + rho * Q1, z), sp.diff(rho**2 * Q2 + rho * Q1, z, 2)),
        "physical_measure": eps * hb * rho,
        "weighted_physical_measure": eps**2 * hb * rho**2,
    }
    p_env = dict(env)
    p_env.update({"fields": fields})
    p_ev = _Eval(p_env, z, qop)
    p_state = _until(pressure, p_ev, p_env, "N")
    if "S" not in p_state or "N" not in p_state:
        raise SourceExtractionError("direct_pressure_functions formula set is incomplete")
    loop = [node for node in pressure.body if isinstance(node, ast.For)]
    if len(loop) != 1 or not isinstance(loop[0].target, ast.Tuple) or [getattr(x, "id", None) for x in loop[0].target.elts] != ["key", "measure", "sign"]:
        raise SourceExtractionError("m13/m15 loop changed")
    loop = loop[0]
    if not isinstance(loop.iter, ast.Tuple) or len(loop.iter.elts) != 2:
        raise SourceExtractionError("m13/m15 entries changed")
    entries = []
    for item in loop.iter.elts:
        if not isinstance(item, ast.Tuple) or len(item.elts) != 3 or not isinstance(item.elts[0], ast.Constant):
            raise SourceExtractionError("m13/m15 entry shape changed")
        entries.append(item)
    if [item.elts[0].value for item in entries] != ["m13", "m15"]:
        raise SourceExtractionError("m13/m15 order changed")
    if _expr(p_ev.run(entries[1].elts[2], p_env)) != -1:
        raise SourceExtractionError("m15 sign changed")
    body_nodes = [node for node in loop.body if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "body" for t in node.targets)]
    if len(body_nodes) != 1:
        raise SourceExtractionError("m15 body assignment changed")
    child = dict(p_env)
    child.update({"key": "m15", "measure": fields["weighted_physical_measure"], "sign": -1, "N": p_state["N"]})
    measure = p_ev.run(entries[1].elts[1], child)
    if sp.simplify(_expr(measure) - fields["weighted_physical_measure"]) != 0:
        raise SourceExtractionError("m15 does not use weighted physical measure")
    m15 = _expr(p_ev.run(_value(body_nodes[0]), child))
    hits = {f"canonical.{name}": _hit(collar, name, state[name]) for name in required}
    hits.update({f"canonical.pressure.{name}": _hit(pressure, name, p_state[name]) for name in ("S", "N")})
    hits["canonical.pressure.m15_body"] = {"function": pressure.name, "assignment": "m15_body", "line": body_nodes[0].lineno, "ast": ast.dump(_value(body_nodes[0]), include_attributes=False), "formula": sp.sstr(m15)}
    return {"collar": state, "pressure": p_state, "m15": m15, "fields": fields}, hits


def _zero(value: sp.Expr, label: str) -> str:
    result = sp.factor(sp.cancel(sp.simplify(value)))
    if result != 0:
        raise SourceExtractionError(f"{label} did not simplify to zero: {result}")
    return "0"


def _ordinary_z1(f: Mapping[str, Any]) -> dict[str, str]:
    z, rho, eps, hb, de, P2 = f["z"], f["rho"], f["eps"], f["hb"], f["de"], f["P2"]
    sigp = f["sigma_prime"]
    R = sp.Function("ratio_of_Z")(z)
    Phi1 = sp.Function("Phi_core_rho_of_Z")(z)
    Psi1 = sp.Function("Psi_core_rho_of_Z")(z)
    cr = (hb - 1) * sigp / (hb * rho)
    L = 1 - de * z**2
    dphi = R * cr * Phi1
    dpsi = R * cr * Psi1
    dq = 2 * z * eps * R * cr * Psi1 / L
    dn = 2 * rho * z * R * cr * Psi1 / L
    A = 2 * rho**2 * z * (1 - hb) * R * Psi1 / L
    dm = -eps**2 * hb * rho**2 * dn
    return {
        "Phi_rhorho": sp.sstr(sp.diff(dphi, z)),
        "Psi_rhorho": sp.sstr(sp.diff(dpsi, z)),
        "Q_rhorho": sp.sstr(sp.diff(dq, z)),
        "N1p_normalized": sp.sstr(sp.diff(dn / P2, z)),
        "N1p_physical": sp.sstr(sp.diff(dn, z)),
        "m15_integrand": sp.sstr(sp.diff(dm, z)),
        "m15_expected_integrand": sp.sstr(sp.diff(eps**2 * A * sigp, z)),
        "all_residuals": _zero(sp.diff(dm - eps**2 * A * sigp, z), "ordinary Z1 m15"),
        "held_constant_in_Z": "epsilon, hb, delta, Pstar2, sigma, sigma_prime, rho, chi_rho",
    }


def prove(forward_source_text, canonical_source_text):
    """Return a JSON-serializable conditional proof report.

    The two positional arguments are the pinned texts for
    :data:`FORWARD_SOURCE_PATH` and :data:`CANONICAL_SOURCE_PATH`.  Any source
    shape change or failed symbolic identity raises :class:`SourceExtractionError`.
    The report is restricted to the common canonical collar interval
    ``[sc_quarter, sc_half]``.
    """

    f, forward_hits, q_template = _forward(forward_source_text)
    c, canonical_hits = _canonical(canonical_source_text, f, q_template)
    z, rho, eps, hb, de, P2 = f["z"], f["rho"], f["eps"], f["hb"], f["de"], f["P2"]
    ratio, phi1, psi1, sigp = f["ratio"], f["phi1"], f["psi1"], f["sigma_prime"]
    L = 1 - de * z**2
    chi_rho = sp.factor(sp.simplify(c["collar"]["Xr"]))
    _zero(chi_rho - (hb - 1) * sigp / (hb * rho), "chi_rho")
    A = 2 * rho**2 * z * (1 - hb) * ratio * psi1 / L
    F, C = f["leading"], c["collar"]
    dphi = _expr(C["phi"][2]) - _expr(F["phi"][2])
    dpsi = _expr(C["psi"][2]) - _expr(F["psi"][2])
    dmean = _expr(C["mean"][2]) - _expr(F["mean"][2])
    dq = _expr(C["Q"][2]) - _expr(F["Q"][2])
    nphys = _expr(c["pressure"]["N"])
    nforward = f["normalized_physical"]
    dn = nphys - nforward
    dm15 = _expr(c["m15"]) - (-eps**2 * hb * rho**2 * nforward)
    residuals = {
        "Phi_rhorho": _zero(dphi - ratio * chi_rho * phi1, "Phi_rhorho"),
        "Psi_rhorho": _zero(dpsi - ratio * chi_rho * psi1, "Psi_rhorho"),
        "mean_rhorho": _zero(dmean, "mean_rhorho"),
        "Q_rhorho": _zero(dq - 2 * z * eps * ratio * chi_rho * psi1 / L, "Q_rhorho"),
        "N1p_physical": _zero(dn - 2 * rho * z * ratio * chi_rho * psi1 / L, "physical N1p"),
        "N1p_normalized": _zero(dn / P2 - 2 * rho * z * ratio * chi_rho * psi1 / (L * P2), "normalized N1p"),
        "m15_ds": _zero(dm15 - eps**2 * A * sigp, "m15 d-sigma source"),
    }
    ordinary = _ordinary_z1(f)
    hits = dict(forward_hits)
    hits.update(canonical_hits)
    return {
        "source_paths": list(SOURCE_PATHS),
        "domain": {"name": "shared_collar_only", "lower": "sc_quarter", "upper": "sc_half"},
        "source_hits": hits,
        "conditions": [
            "shared real leading Phi, Psi, and mean fields on the common collar domain",
            "canonical ratio/Psi/mean history is represented by the extracted AST bodies",
            "canonical own.radial_Q uses the radial operator template extracted from forward Q; this is an explicit conditional operator contract, not a cross-owner proof",
            "rho = 4 exp(hb s), d rho/ds = hb rho, chi_rho = (hb - 1) sigma_prime/(hb rho)",
            "epsilon, hb, delta, Pstar2, sigma, sigma_prime, rho, and chi_rho are held constant under ordinary Z1",
        ],
        "identities": {
            "canonical_minus_forward_Phi_rhorho": {"residual": residuals["Phi_rhorho"], "difference": sp.sstr(ratio * chi_rho * phi1)},
            "canonical_minus_forward_Psi_rhorho": {"residual": residuals["Psi_rhorho"], "difference": sp.sstr(ratio * chi_rho * psi1)},
            "canonical_minus_forward_mean_rhorho": {"residual": residuals["mean_rhorho"], "difference": "0"},
            "canonical_minus_forward_Q_rhorho": {"residual": residuals["Q_rhorho"], "difference": sp.sstr(2 * z * eps * ratio * chi_rho * psi1 / L)},
            "canonical_minus_forward_N1p_physical": {"residual": residuals["N1p_physical"], "difference": sp.sstr(2 * rho * z * ratio * chi_rho * psi1 / L)},
            "canonical_minus_forward_N1p_normalized": {"residual": residuals["N1p_normalized"], "difference": sp.sstr(2 * rho * z * ratio * chi_rho * psi1 / (L * P2))},
            "canonical_minus_forward_m15_ds": {"residual": residuals["m15_ds"], "difference": sp.sstr(eps**2 * A * sigp)},
        },
        "ordinary_Z1": ordinary,
        "scope_exclusions": [
            "no cross-owner same-source identity",
            "no numerical value, graph evaluation, or saved-data read",
            "no claim that forward six finite terms equal canonical inner plus lower",
            "no error, tail, or cancellation proof",
            "no assertion on the non-common q = 0..2 range",
        ],
    }


__all__ = ["FORWARD_SOURCE_PATH", "CANONICAL_SOURCE_PATH", "SOURCE_PATHS", "SourceExtractionError", "prove"]

