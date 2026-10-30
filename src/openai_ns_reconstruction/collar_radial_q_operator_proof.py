"""Strict source-text proof of the original collar radial ``Q`` operator.

``prove`` accepts only the four requested source strings.  It performs no
imports, filesystem reads, owner construction, numerical evaluation, or test
execution.  AST gates reject source changes before the symbolic C2 comparison.
The result is an operator identity only: source-history and cross-owner
identity remain explicitly unproved, and no finite approximation is promoted
to a function.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
import textwrap
from typing import Any, Iterable

import sympy as sp


class CanonicalSourceError(ValueError):
    """An input source text changed outside the admitted canonical AST."""


RADIAL_BODY = """
    a=c2.C2Algebra(g);z=c2.C2Function(g.symbol('Z'),g.one,g.zero)
    one=a.fixed(1);d=a.add(one,a.neg(a.mul(z,z)))
    de=a.fixed(delta);L=a.add(one,a.neg(a.mul(de,a.mul(z,z))))
    u=c2.C2Function(U.value,U.Z,U.ZZ);mean=c2.C2Function(m.value,m.Z,m.ZZ)
    mz=c2.C2Function(m.Z,m.ZZ,m.ZZZ)
    numerator=a.add(a.scale(a.mul(z,u),2),
        a.neg(a.mul(a.add(one,a.neg(de)),a.mul(z,mean))),a.neg(a.mul(d,mz)))
    return a.div(numerator,L)
"""
FORWARD_DL = """
    d,L=1-z*z,1-z*z*de
"""
FORWARD_Q = """
    Q=[(u*(2*z)-m*(z*(1-de))-m.dz()*d)/L for u,m in zip(U,M)]
"""
C2_METHODS = {
    "fixed": """
        q=self.g.constant(value) if isinstance(value,(str,int,Fraction)) else value
        return C2Function(q,self.g.zero,self.g.zero)
    """,
    "add": """
        return C2Function(*(self.g.add(*(getattr(q,key) for q in rows)) for key in ('value','Z','ZZ')))
    """,
    "neg": """
        return C2Function(*(self.g.neg(getattr(q,key)) for key in ('value','Z','ZZ')))
    """,
    "scale": """
        return self.mul(q,self.fixed(v))
    """,
    "mul": """
        g=self.g
        return C2Function(g.mul(a.value,b.value),g.add(g.mul(a.Z,b.value),g.mul(a.value,b.Z)),
            g.add(g.mul(a.ZZ,b.value),g.mul(g.constant(2),a.Z,b.Z),g.mul(a.value,b.ZZ)))
    """,
    "div": """
        g=self.g;divide=lambda a:g.quotient(a,d.value,'same actual positive source denominator')
        v=divide(n.value);z=divide(g.sub(n.Z,g.mul(v,d.Z)))
        zz=divide(g.sub(n.ZZ,g.add(g.mul(g.constant(2),z,d.Z),g.mul(v,d.ZZ))))
        return C2Function(v,z,zz)
    """,
}
LAURENT_DZ = """
    return Laurent({p:type(v)(v.ctx,[(k+1)*v[k+1] for k in range(v.order)]) for p,v in self.rows.items()})
"""


def _dump(node: ast.AST) -> str:
    return ast.dump(node, annotate_fields=True, include_attributes=False)


def _parse(text: str, label: str) -> ast.Module:
    if not isinstance(text, str):
        raise TypeError(f"{label} must be source text")
    try:
        return ast.parse(text, filename=label, mode="exec")
    except SyntaxError as exc:
        raise CanonicalSourceError(f"invalid {label}: {exc}") from exc


def _body_without_docstring(body: list[ast.stmt]) -> list[ast.stmt]:
    if body and isinstance(body[0], ast.Expr):
        value = body[0].value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            return body[1:]
    return body


def _fragment(source: str, label: str) -> list[ast.stmt]:
    body = textwrap.indent(textwrap.dedent(source).strip("\n"), "    ")
    module = ast.parse("def _fragment():\n" + body, filename=label, mode="exec")
    function = module.body[0]
    assert isinstance(function, ast.FunctionDef)
    return _body_without_docstring(function.body)


def _function(module: ast.Module, name: str, label: str) -> ast.FunctionDef:
    matches = [node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == name]
    if len(matches) != 1:
        raise CanonicalSourceError(f"{label}: expected one top-level {name}")
    return matches[0]


def _class(module: ast.Module, name: str, label: str) -> ast.ClassDef:
    matches = [node for node in module.body if isinstance(node, ast.ClassDef) and node.name == name]
    if len(matches) != 1:
        raise CanonicalSourceError(f"{label}: expected one top-level {name}")
    return matches[0]


def _plain_args(function: ast.FunctionDef, names: tuple[str, ...], label: str) -> None:
    args = function.args
    if args.posonlyargs or args.vararg or args.kwonlyargs or args.kwarg or args.defaults or args.kw_defaults:
        raise CanonicalSourceError(f"{label}: argument form changed")
    if tuple(arg.arg for arg in args.args) != names:
        raise CanonicalSourceError(f"{label}: argument names changed")


def _method(class_node: ast.ClassDef, name: str, label: str) -> ast.FunctionDef:
    matches = [node for node in class_node.body if isinstance(node, ast.FunctionDef) and node.name == name]
    if len(matches) != 1:
        raise CanonicalSourceError(f"{label}: expected one {name} method")
    return matches[0]


def _method_args(function: ast.FunctionDef, positional: tuple[str, ...], vararg: str | None, label: str) -> None:
    args = function.args
    actual_vararg = args.vararg.arg if args.vararg else None
    if args.posonlyargs or args.kwonlyargs or args.kwarg or args.defaults or args.kw_defaults:
        raise CanonicalSourceError(f"{label}: argument form changed")
    if tuple(arg.arg for arg in args.args) != positional or actual_vararg != vararg:
        raise CanonicalSourceError(f"{label}: argument names changed")


def _validate_radial(own: ast.Module) -> ast.FunctionDef:
    function = _function(own, "radial_Q", "own_text")
    _plain_args(function, ("g", "delta", "U", "m"), "radial_Q")
    if [_dump(node) for node in _body_without_docstring(function.body)] != [_dump(node) for node in _fragment(RADIAL_BODY, "canonical radial_Q")]:
        raise CanonicalSourceError("radial_Q: canonical assignments or return changed")
    return function


def _validate_forward(forward: ast.Module) -> tuple[ast.Assign, ast.Assign]:
    function = _function(forward, "leading_regular_frame", "forward_text")
    _plain_args(function, ("values", "comparison", "z", "rho", "x", "de", "eps", "P2", "sigma", "hb", "axis_j"), "leading_regular_frame")
    body = _body_without_docstring(function.body)
    dl = [node for node in body if isinstance(node, ast.Assign) and len(node.targets) == 1
          and isinstance(node.targets[0], (ast.Tuple, ast.List))
          and [item.id if isinstance(item, ast.Name) else None for item in node.targets[0].elts] == ["d", "L"]]
    q = [node for node in body if isinstance(node, ast.Assign) and len(node.targets) == 1
         and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "Q"]
    if len(dl) != 1 or len(q) != 1:
        raise CanonicalSourceError("leading_regular_frame: d/L or Q assignment changed")
    if _dump(dl[0]) != _dump(_fragment(FORWARD_DL, "canonical forward d/L")[0]):
        raise CanonicalSourceError("leading_regular_frame: d/L changed")
    if _dump(q[0]) != _dump(_fragment(FORWARD_Q, "canonical forward Q")[0]):
        raise CanonicalSourceError("leading_regular_frame: Q changed")
    if body.index(q[0]) != body.index(dl[0]) + 1:
        raise CanonicalSourceError("leading_regular_frame: Q must follow d/L")
    return dl[0], q[0]


def _validate_c2(c2: ast.Module) -> None:
    function = _class(c2, "C2Function", "c2_text")
    fields = [node.target.id for node in function.body if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)]
    if fields != ["value", "Z", "ZZ"]:
        raise CanonicalSourceError("C2Function shape changed")
    algebra = _class(c2, "C2Algebra", "c2_text")
    signatures = {
        "fixed": (("self", "value"), None), "add": (("self",), "rows"),
        "neg": (("self", "q"), None), "scale": (("self", "q", "v"), None),
        "mul": (("self", "a", "b"), None), "div": (("self", "n", "d"), None),
    }
    for name, source in C2_METHODS.items():
        method = _method(algebra, name, "C2Algebra")
        positional, vararg = signatures[name]
        _method_args(method, positional, vararg, f"C2Algebra.{name}")
        if [_dump(node) for node in _body_without_docstring(method.body)] != [_dump(node) for node in _fragment(source, f"canonical C2Algebra.{name}")]:
            raise CanonicalSourceError(f"C2Algebra.{name}: ordinary derivative formula changed")


def _validate_laurent(laurent: ast.Module) -> None:
    method = _method(_class(laurent, "Laurent", "laurent_text"), "dz", "Laurent")
    _plain_args(method, ("self",), "Laurent.dz")
    if [_dump(node) for node in _body_without_docstring(method.body)] != [_dump(node) for node in _fragment(LAURENT_DZ, "canonical Laurent.dz")]:
        raise CanonicalSourceError("Laurent.dz: Taylor derivative map changed")


@dataclass(frozen=True)
class _C2Jet:
    value: sp.Expr
    Z: sp.Expr
    ZZ: sp.Expr


class _C2:
    @staticmethod
    def fixed(value: sp.Expr) -> _C2Jet:
        return _C2Jet(value, sp.Integer(0), sp.Integer(0))

    @staticmethod
    def add(*rows: _C2Jet) -> _C2Jet:
        return _C2Jet(*(sp.Add(*(getattr(row, key) for row in rows)) for key in ("value", "Z", "ZZ")))

    @staticmethod
    def neg(row: _C2Jet) -> _C2Jet:
        return _C2Jet(-row.value, -row.Z, -row.ZZ)

    @staticmethod
    def scale(row: _C2Jet, value: sp.Expr) -> _C2Jet:
        return _C2.mul(row, _C2.fixed(value))

    @staticmethod
    def mul(left: _C2Jet, right: _C2Jet) -> _C2Jet:
        return _C2Jet(left.value * right.value,
                       left.Z * right.value + left.value * right.Z,
                       left.ZZ * right.value + 2 * left.Z * right.Z + left.value * right.ZZ)

    @staticmethod
    def div(numerator: _C2Jet, denominator: _C2Jet) -> _C2Jet:
        value = numerator.value / denominator.value
        first = (numerator.Z - value * denominator.Z) / denominator.value
        second = (numerator.ZZ - 2 * first * denominator.Z - value * denominator.ZZ) / denominator.value
        return _C2Jet(value, first, second)


def _validated_current_q() -> tuple[_C2Jet, sp.Symbol, sp.Symbol]:
    z = sp.Symbol("z", real=True)
    delta = sp.Symbol("delta", real=True)
    zj = _C2Jet(z, sp.Integer(1), sp.Integer(0))
    one = _C2.fixed(sp.Integer(1))
    d = _C2.add(one, _C2.neg(_C2.mul(zj, zj)))
    L = _C2.add(one, _C2.neg(_C2.mul(_C2.fixed(delta), _C2.mul(zj, zj))))
    u = _C2Jet(sp.Symbol("u0"), sp.Symbol("u1"), sp.Symbol("u2"))
    m = _C2Jet(sp.Symbol("m0"), sp.Symbol("m1"), sp.Symbol("m2"))
    mz = _C2Jet(sp.Symbol("m1"), sp.Symbol("m2"), sp.Symbol("m3"))
    numerator = _C2.add(_C2.scale(_C2.mul(zj, u), sp.Integer(2)),
                         _C2.neg(_C2.mul(_C2.add(one, _C2.neg(_C2.fixed(delta))), _C2.mul(zj, m))),
                         _C2.neg(_C2.mul(d, mz)))
    return _C2.div(numerator, L), z, delta


def _forward_rows() -> tuple[tuple[sp.Expr, sp.Expr, sp.Expr], sp.Symbol, sp.Symbol]:
    z = sp.Symbol("z", real=True)
    delta = sp.Symbol("delta", real=True)
    u = sp.Function("u")(z)
    m = sp.Function("m")(z)
    L = 1 - delta * z**2
    q = (2 * z * u - (1 - delta) * z * m - (1 - z**2) * sp.diff(m, z)) / L
    qz, qzz = sp.diff(q, z), sp.diff(q, z, 2)
    u0, u1, u2 = sp.symbols("u0 u1 u2")
    m0, m1, m2, m3 = sp.symbols("m0 m1 m2 m3")
    replacement = {
        u: u0, sp.diff(u, z): u1, sp.diff(u, z, 2): u2,
        m: m0, sp.diff(m, z): m1, sp.diff(m, z, 2): m2, sp.diff(m, z, 3): m3,
    }
    return tuple(sp.cancel(row.xreplace(replacement)) for row in (q, qz, qzz)), z, delta


def _difference(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.cancel(sp.together(left - right)))


def _text(value: sp.Expr) -> str:
    return sp.sstr(sp.factor(value))


def prove(own_text: str, forward_text: str, c2_text: str, laurent_text: str) -> dict[str, Any]:
    """Strictly validate four source texts and prove the C2 Q rows agree."""
    own = _parse(own_text, "own_text")
    forward = _parse(forward_text, "forward_text")
    c2 = _parse(c2_text, "c2_text")
    laurent = _parse(laurent_text, "laurent_text")
    _validate_radial(own)
    _validate_forward(forward)
    _validate_c2(c2)
    _validate_laurent(laurent)

    current, z, delta = _validated_current_q()
    forward_rows, forward_z, forward_delta = _forward_rows()
    if (z, delta) != (forward_z, forward_delta):
        raise CanonicalSourceError("symbolic variable alignment changed")
    differences = tuple(_difference(left, right) for left, right in zip((current.value, current.Z, current.ZZ), forward_rows))
    if any(value != 0 for value in differences):
        raise CanonicalSourceError("validated Q operators differ")
    L = 1 - delta * z**2
    return {
        "passed": True,
        "ast_validation": {
            "canonical_radial_Q": True,
            "forward_d_and_L": True,
            "forward_Q": True,
            "C2Function_shape": "(value, Z, ZZ)",
            "C2Algebra_fixed_add_neg_scale_mul_div": True,
            "Laurent_dz": True,
        },
        "owner": {"source_function": "radial_Q", "placeholder_owner_used": False, "v6_placeholder_used": False},
        "current_C2_rows": [_text(row) for row in (current.value, current.Z, current.ZZ)],
        "forward_ordinary_Z_rows": [_text(row) for row in forward_rows],
        "differences_by_order": ["0", "0", "0"],
        "difference_sympy": [_text(value) for value in differences],
        "L": {
            "expression": _text(L),
            "required": "1 - delta*z**2 != 0",
            "sufficient_real_condition": "abs(z) < 1 and 0 <= delta < 1",
            "canonical_stronger_condition": "abs(z) < 1 and 0 < delta < 1/2 implies L > 0",
        },
        "Taylor_to_ordinary_Z": {
            "coefficient_k_to_derivative": "f^(k) = k! * coefficient[k]",
            "Laurent_dz_row": "[(k+1)*v[k+1] for k in range(v.order)]",
            "C3_rows": "(value, Z, ZZ, ZZZ) are ordinary derivatives",
        },
        "limitations": {
            "source_history_identity_proved": False,
            "cross_owner_identity_proved": False,
            "finite_approximation_identified_as_function": False,
            "numerical_values_evaluated": False,
        },
    }


__all__ = ["CanonicalSourceError", "prove"]
