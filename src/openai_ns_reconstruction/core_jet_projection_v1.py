"""Static source proof for the exact first-order core-jet projection.

The public entry point :func:`build` accepts the complete source strings of
the high/source-matrix and low/matrix leaves.  It only parses those strings
with :mod:`ast`; it never imports, opens, evaluates, or executes a producer.
The source anchors are deliberately narrow: they require the factorial source
jet, convolution product, reciprocal recurrence, derivative, relative-F0
recurrence, low ``Dual`` product/reciprocal, and the explicit low ``fz``
second-slot formula.

After the anchors pass, a tiny exact rational-function algebra proves the
projection theorem in ``Q(a,b,...)[h]/(h**2)``.  This is a scalar,
commutative-ring identity conditional on nonzero denominators and on the two
leaves receiving the same source atoms.  It does not claim equality of
interval endpoints, physical-field ownership, a sum over fields, or a
numerical evaluation.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Iterable


SCHEMA = "core_jet_projection_v1"
HIGH_LABEL = "n1_core_low_source_matrix_v2.py"
LOW_LABEL = "n1_core_low_matrix_v3.py"

ASSUMPTIONS = (
    "The two inputs are complete source strings supplied by the caller.",
    "Both strings are inspected only with ast.parse and AST shape predicates.",
    "The scalar coefficients lie in a commutative ring; displayed quotients use nonzero denominators.",
    "The high coefficients are factorial-normalised Taylor coefficients.",
    "The low Dual derivative slot is the ordinary first derivative with respect to the same h/Z.",
    "High and low source rows denote the same input atoms before projection.",
)

UNPROVED = (
    "interval endpoints are not asserted equal",
    "no physical field is materialised or owned",
    "no sum over fields or matrix rows is claimed",
    "no numerical interval enclosure or endpoint equality is computed",
    "source provenance equality is assumed, not inferred from names",
    "the source AST anchors are partial shape checks, not a full verification of every constructor, operator or ell definition",
)

SCOPE_NOTES = (
    "The low-matrix prose/docstring description of Phi/F0 differentiation is semantically misleading; this receipt trusts the explicit fz AST formula, not that prose.",
    "The proved fz identity is for F=Phi*C with C0=1 and C1=ell0; Phi_Z+ell*Phi is the normalized derivative of F0*Phi, not the derivative of Phi/F0.",
)

HARD_TO_CERTIFY = (
    "AST shape cannot prove that two independently supplied source atoms have the same physical owner.",
    "AST shape cannot prove interval endpoint equality or directed-rounding equality.",
    "The scalar proof abstracts context arithmetic to a commutative exact ring and therefore does not certify interval implementation bounds.",
)


class SourceIdentityError(ValueError):
    """Reserved for callers that want to classify malformed source text."""


class _AnchorError(ValueError):
    pass


# ---------------------------------------------------------------------------
# A small exact scalar algebra: Q(a,b,...), represented as rational functions.


Monomial = tuple[str, ...]


def _sorted_terms(values: dict[Monomial, Fraction]) -> tuple[tuple[Monomial, Fraction], ...]:
    return tuple(
        sorted(
            ((monomial, coefficient) for monomial, coefficient in values.items() if coefficient),
            key=lambda item: (len(item[0]), item[0]),
        )
    )


@dataclass(frozen=True)
class _Poly:
    terms: tuple[tuple[Monomial, Fraction], ...] = ()

    @staticmethod
    def zero() -> "_Poly":
        return _Poly(())

    @staticmethod
    def one() -> "_Poly":
        return _Poly((((), Fraction(1)),))

    @staticmethod
    def atom(name: str) -> "_Poly":
        return _Poly((( (name,), Fraction(1)),))

    @staticmethod
    def constant(value: int | Fraction) -> "_Poly":
        coefficient = Fraction(value)
        return _Poly(()) if coefficient == 0 else _Poly((((), coefficient),))

    def _dict(self) -> dict[Monomial, Fraction]:
        return dict(self.terms)

    def __add__(self, other: "_Poly") -> "_Poly":
        values = self._dict()
        for monomial, coefficient in other.terms:
            values[monomial] = values.get(monomial, Fraction(0)) + coefficient
        return _Poly(_sorted_terms(values))

    def __neg__(self) -> "_Poly":
        return _Poly(tuple((m, -c) for m, c in self.terms))

    def __sub__(self, other: "_Poly") -> "_Poly":
        return self + (-other)

    def __mul__(self, other: "_Poly") -> "_Poly":
        values: dict[Monomial, Fraction] = {}
        for left_monomial, left_coefficient in self.terms:
            for right_monomial, right_coefficient in other.terms:
                monomial = tuple(sorted(left_monomial + right_monomial))
                values[monomial] = values.get(monomial, Fraction(0)) + left_coefficient * right_coefficient
        return _Poly(_sorted_terms(values))

    def __pow__(self, power: int) -> "_Poly":
        if power < 0:
            raise ValueError("negative polynomial powers belong to the rational layer")
        result = _Poly.one()
        for _ in range(power):
            result = result * self
        return result

    @property
    def is_zero(self) -> bool:
        return not self.terms

    def text(self) -> str:
        if not self.terms:
            return "0"
        pieces: list[str] = []
        for monomial, coefficient in self.terms:
            sign = "-" if coefficient < 0 else "+"
            magnitude = -coefficient if coefficient < 0 else coefficient
            if magnitude.denominator == 1:
                number = str(magnitude.numerator)
            else:
                number = f"{magnitude.numerator}/{magnitude.denominator}"
            factors: list[str] = []
            for name in monomial:
                if factors and factors[-1].startswith(name + "^"):
                    continue
                count = monomial.count(name)
                factors.append(name if count == 1 else f"{name}^{count}")
            if factors and magnitude == 1:
                body = "*".join(factors)
            elif factors:
                body = number + "*" + "*".join(factors)
            else:
                body = number
            if not pieces:
                pieces.append(("-" if sign == "-" else "") + body)
            else:
                pieces.append((" - " if sign == "-" else " + ") + body)
        return "".join(pieces)


def _poly(value: _Poly | int | Fraction) -> _Poly:
    return value if isinstance(value, _Poly) else _Poly.constant(value)


@dataclass(frozen=True)
class _Rat:
    numerator: _Poly
    denominator: _Poly

    @staticmethod
    def constant(value: int | Fraction) -> "_Rat":
        return _Rat(_Poly.constant(value), _Poly.one())

    @staticmethod
    def atom(name: str) -> "_Rat":
        return _Rat(_Poly.atom(name), _Poly.one())

    def __add__(self, other: "_Rat") -> "_Rat":
        return _Rat(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator,
        )

    def __neg__(self) -> "_Rat":
        return _Rat(-self.numerator, self.denominator)

    def __sub__(self, other: "_Rat") -> "_Rat":
        return self + (-other)

    def __mul__(self, other: "_Rat") -> "_Rat":
        return _Rat(self.numerator * other.numerator, self.denominator * other.denominator)

    def __rmul__(self, other: int | Fraction) -> "_Rat":
        return _Rat(_poly(other) * self.numerator, self.denominator)

    def reciprocal(self) -> "_Rat":
        if self.numerator.is_zero:
            raise ZeroDivisionError("formal scalar denominator is zero")
        return _Rat(self.denominator, self.numerator)

    def __truediv__(self, other: "_Rat") -> "_Rat":
        return self * other.reciprocal()

    def __pow__(self, power: int) -> "_Rat":
        if power < 0:
            return self.reciprocal() ** (-power)
        return _Rat(self.numerator ** power, self.denominator ** power)

    def equal(self, other: "_Rat") -> bool:
        return self.numerator * other.denominator == other.numerator * self.denominator

    def text(self) -> str:
        numerator = self.numerator.text()
        denominator = self.denominator.text()
        if denominator == "1":
            return numerator
        return f"({numerator})/({denominator})"


def _r(value: int | Fraction) -> _Rat:
    return _Rat.constant(value)


def _v(name: str) -> _Rat:
    return _Rat.atom(name)


@dataclass(frozen=True)
class _Jet:
    """Factorial-normalised truncated Taylor coefficients."""

    coefficients: tuple[_Rat, ...]

    def add(self, other: "_Jet") -> "_Jet":
        order = min(len(self.coefficients), len(other.coefficients))
        return _Jet(tuple(self.coefficients[i] + other.coefficients[i] for i in range(order)))

    def neg(self) -> "_Jet":
        return _Jet(tuple(-value for value in self.coefficients))

    def mul(self, other: "_Jet") -> "_Jet":
        order = min(len(self.coefficients), len(other.coefficients))
        values: list[_Rat] = []
        for index in range(order):
            value = _r(0)
            for left in range(index + 1):
                value = value + self.coefficients[left] * other.coefficients[index - left]
            values.append(value)
        return _Jet(tuple(values))

    def square(self) -> "_Jet":
        return self.mul(self)

    def reciprocal(self) -> "_Jet":
        if not self.coefficients:
            raise ZeroDivisionError("empty formal jet")
        values = [_r(1) / self.coefficients[0]]
        for index in range(1, len(self.coefficients)):
            value = _r(0)
            for left in range(1, index + 1):
                value = value + self.coefficients[left] * values[index - left]
            values.append(-value / self.coefficients[0])
        return _Jet(tuple(values))

    def div(self, other: "_Jet") -> "_Jet":
        return self.mul(other.reciprocal())

    def differentiate(self) -> "_Jet":
        if len(self.coefficients) < 2:
            raise ValueError("a derivative needs a first Taylor coefficient")
        return _Jet(tuple((index + 1) * value for index, value in enumerate(self.coefficients[1:])))

    def project(self) -> tuple[_Rat, _Rat]:
        if len(self.coefficients) < 2:
            raise ValueError("projection needs order one")
        return self.coefficients[0], self.coefficients[1]


@dataclass(frozen=True)
class _DualPair:
    value: _Rat
    derivative: _Rat

    def add(self, other: "_DualPair") -> "_DualPair":
        return _DualPair(self.value + other.value, self.derivative + other.derivative)

    def neg(self) -> "_DualPair":
        return _DualPair(-self.value, -self.derivative)

    def mul(self, other: "_DualPair") -> "_DualPair":
        return _DualPair(
            self.value * other.value,
            self.derivative * other.value + self.value * other.derivative,
        )

    def square(self) -> "_DualPair":
        return self.mul(self)

    def reciprocal(self) -> "_DualPair":
        return _DualPair(
            self.value.reciprocal(),
            -self.derivative / (self.value ** 2),
        )

    def div(self, other: "_DualPair") -> "_DualPair":
        return self.mul(other.reciprocal())

    def pair(self) -> tuple[_Rat, _Rat]:
        return self.value, self.derivative


def _pair_text(pair: tuple[_Rat, _Rat]) -> tuple[str, str]:
    return tuple(value.text() for value in pair)  # type: ignore[return-value]


def _same_pair(left: tuple[_Rat, _Rat], right: tuple[_Rat, _Rat]) -> bool:
    return all(a.equal(b) for a, b in zip(left, right))


def _exact_algebra_report() -> dict[str, Any]:
    a = _Jet((_v("a0"), _v("a1")))
    b = _Jet((_v("b0"), _v("b1")))
    da = _DualPair(_v("a0"), _v("a1"))
    db = _DualPair(_v("b0"), _v("b1"))
    operations: dict[str, Any] = {}

    cases = (
        ("add", a.add(b).project(), da.add(db).pair(), "(a0+b0, a1+b1)"),
        ("neg", a.neg().project(), da.neg().pair(), "(-a0, -a1)"),
        ("mul", a.mul(b).project(), da.mul(db).pair(), "(a0*b0, a0*b1+a1*b0)"),
        ("square", a.square().project(), da.square().pair(), "(a0^2, 2*a0*a1)"),
        ("reciprocal", a.reciprocal().project(), da.reciprocal().pair(), "(a0^-1, -a1*a0^-2)"),
        (
            "div",
            a.div(b).project(),
            da.div(db).pair(),
            "(a0/b0, (a1*b0-a0*b1)/b0^2)",
        ),
    )
    for name, high_pair, dual_pair, closed_form in cases:
        operations[name] = {
            "ring_formula": closed_form,
            "high_projection": _pair_text(high_pair),
            "dual_pair": _pair_text(dual_pair),
            "matched": _same_pair(high_pair, dual_pair),
        }

    second = _Jet((_v("a0"), _v("a1"), _v("a2")))
    high_derivative = second.differentiate().project()
    dual_derivative = (_v("a1"), _r(2) * _v("a2"))
    derivative = {
        "ring_formula": "pi(d/dh(a0+a1*h+a2*h^2)) = (a1, 2*a2)",
        "high_projection": _pair_text(high_derivative),
        "dual_pair": _pair_text(dual_derivative),
        "matched": _same_pair(high_derivative, dual_derivative),
    }

    ell0, ell1 = _v("ell0"), _v("ell1")
    c0 = _r(1)
    c1 = ell0 * c0
    c2 = (ell0 * c1 + ell1 * c0) / _r(2)
    phi0, phi1, phi2 = _v("Phi0"), _v("Phi1"), _v("Phi2")
    phi = _Jet((phi0, phi1, phi2 / _r(2)))
    relative = _Jet((c0, c1, c2))
    product = phi.mul(relative)
    F0, F1, F2_coeff = product.coefficients
    F2_ordinary = _r(2) * F2_coeff
    recurrence = {
        "equation": "C'(h)=ell(h)*C(h), C0=1",
        "C0": c0.text(),
        "C1": c1.text(),
        "two_C2": (ell1 + ell0 ** 2).text(),
        "C2": c2.text(),
        "matched": c0.equal(_r(1)) and c1.equal(ell0) and c2.equal((ell1 + ell0 ** 2) / _r(2)),
    }
    fz_dual = _DualPair(F1, F2_ordinary)
    expected_fz = _DualPair(
        phi1 + ell0 * phi0,
        phi2 + _r(2) * ell0 * phi1 + (ell1 + ell0 ** 2) * phi0,
    )
    f_product = {
        "F0": F0.text(),
        "F1": F1.text(),
        "F2_coefficient": F2_coeff.text(),
        "F2_ordinary": F2_ordinary.text(),
        "expected_F0": "Phi0",
        "expected_F1": "Phi1+ell0*Phi0",
        "expected_F2_ordinary": "Phi2+2*ell0*Phi1+(ell1+ell0^2)*Phi0",
        "matched": F0.equal(phi0)
        and F1.equal(phi1 + ell0 * phi0)
        and F2_ordinary.equal(phi2 + _r(2) * ell0 * phi1 + (ell1 + ell0 ** 2) * phi0),
        "fz_dual": {
            "computed": _pair_text(fz_dual.pair()),
            "expected": _pair_text(expected_fz.pair()),
            "matched": _same_pair(fz_dual.pair(), expected_fz.pair()),
        },
    }
    all_matched = all(item["matched"] for item in operations.values()) and derivative["matched"] and recurrence["matched"] and f_product["matched"] and f_product["fz_dual"]["matched"]
    return {
        "ring": "Q(a0,a1,b0,b1,ell0,ell1,Phi0,Phi1,Phi2)[h]/(h^2)",
        "projection": "pi(a0+a1*h)=(a0,a1)",
        "operations": operations,
        "differentiate": derivative,
        "relative_F0_recurrence": recurrence,
        "F_product": f_product,
        "exact": all_matched,
    }


# ---------------------------------------------------------------------------
# AST extraction and narrow implementation anchors.


def _ast_dump(node: ast.AST | None) -> str | None:
    return None if node is None else ast.dump(node, annotate_fields=True, include_attributes=False)


def _dotted(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        prefix = _dotted(node.value)
        return None if prefix is None else prefix + "." + node.attr
    return None


def _name(node: ast.AST | None, expected: str) -> bool:
    return isinstance(node, ast.Name) and node.id == expected


def _integer(node: ast.AST | None, expected: int) -> bool:
    return isinstance(node, ast.Constant) and type(node.value) is int and node.value == expected


def _subscript(node: ast.AST | None, root: str, index: ast.AST | None = None) -> bool:
    if not isinstance(node, ast.Subscript) or not _name(node.value, root):
        return False
    slice_node = node.slice
    if hasattr(ast, "Index") and isinstance(slice_node, ast.Index):
        slice_node = slice_node.value
    return index is None or _shape(slice_node) == _shape(index)


def _shape(node: ast.AST | None) -> tuple[Any, ...] | None:
    if node is None:
        return None
    if isinstance(node, ast.Name):
        return ("name", node.id)
    if isinstance(node, ast.Constant) and type(node.value) in (int, float, str):
        return ("constant", node.value)
    if isinstance(node, ast.Attribute):
        return ("attribute", _dotted(node))
    if isinstance(node, ast.Subscript):
        slice_node = node.slice
        if hasattr(ast, "Index") and isinstance(slice_node, ast.Index):
            slice_node = slice_node.value
        return ("subscript", _shape(node.value), _shape(slice_node))
    if isinstance(node, ast.BinOp):
        operators = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/", ast.Pow: "**"}
        op = next((symbol for kind, symbol in operators.items() if isinstance(node.op, kind)), None)
        return None if op is None else ("bin", op, _shape(node.left), _shape(node.right))
    if isinstance(node, ast.UnaryOp):
        op = "-" if isinstance(node.op, ast.USub) else "+" if isinstance(node.op, ast.UAdd) else None
        return None if op is None else ("unary", op, _shape(node.operand))
    if isinstance(node, ast.Call):
        return ("call", _dotted(node.func), tuple(_shape(arg) for arg in node.args))
    return (type(node).__name__,)


def _args(node: ast.FunctionDef | ast.AsyncFunctionDef) -> tuple[str, ...]:
    positional = list(getattr(node.args, "posonlyargs", ())) + list(node.args.args)
    return tuple(item.arg for item in positional)


def _signature(node: ast.AST | None, expected: tuple[str, ...]) -> bool:
    if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        return False
    return (
        _args(node) == expected
        and not node.args.vararg
        and not node.args.kwarg
        and not node.args.kwonlyargs
        and not node.args.defaults
        and not node.args.kw_defaults
    )


def _functions(tree: ast.AST | None, name: str) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    if tree is None:
        return []
    return [node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]


def _class(tree: ast.AST | None, name: str) -> ast.ClassDef | None:
    if tree is None:
        return None
    classes = [node for node in ast.walk(tree) if isinstance(node, ast.ClassDef) and node.name == name]
    return classes[0] if len(classes) == 1 else None


def _method(tree: ast.AST | None, class_name: str, method_name: str) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    owner = _class(tree, class_name)
    if owner is None:
        return None
    methods = [node for node in owner.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == method_name]
    return methods[0] if len(methods) == 1 else None


def _single_function(tree: ast.AST | None, name: str) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    found = _functions(tree, name)
    return found[0] if len(found) == 1 else None


def _assignments(scope: ast.AST | None, target: str) -> list[ast.Assign]:
    if scope is None:
        return []
    return [
        node
        for node in ast.walk(scope)
        if isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and _name(node.targets[0], target)
    ]


def _calls(node: ast.AST | None, function: str) -> list[ast.Call]:
    if node is None:
        return []
    return [item for item in ast.walk(node) if isinstance(item, ast.Call) and _dotted(item.func) == function]


def _contains(node: ast.AST | None, predicate: Any) -> bool:
    return any(predicate(item) for item in ast.walk(node)) if node is not None else False


def _is_add_one(node: ast.AST | None, variable: str) -> bool:
    return isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add) and ((_name(node.left, variable) and _integer(node.right, 1)) or (_integer(node.left, 1) and _name(node.right, variable)))


def _is_range_plus_one(node: ast.AST | None, variable: str) -> bool:
    return isinstance(node, ast.Call) and _dotted(node.func) == "range" and len(node.args) == 1 and _is_add_one(node.args[0], variable)


def _is_range_from_one(node: ast.AST | None, variable: str) -> bool:
    return isinstance(node, ast.Call) and _dotted(node.func) == "range" and len(node.args) == 2 and _integer(node.args[0], 1) and _is_add_one(node.args[1], variable)


def _is_factorial_call(node: ast.AST | None, variable: str) -> bool:
    return isinstance(node, ast.Call) and _dotted(node.func) in {"math.factorial", "factorial"} and len(node.args) == 1 and _name(node.args[0], variable)


def _is_n_denominator(node: ast.AST | None) -> bool:
    return _name(node, "n") or (
        isinstance(node, ast.Call)
        and _dotted(node.func) == "context.mpf"
        and len(node.args) == 1
        and _name(node.args[0], "n")
    )


def _is_index_difference(node: ast.AST | None, outer: str, inner: str) -> bool:
    return isinstance(node, ast.BinOp) and isinstance(node.op, ast.Sub) and _name(node.left, outer) and _name(node.right, inner)


def _is_convolution_term(node: ast.AST | None, left_root: str, left_index: str, right_root: str, right_outer: str, right_inner: str) -> bool:
    if not isinstance(node, ast.BinOp) or not isinstance(node.op, ast.Mult):
        return False
    left, right = node.left, node.right
    return (_subscript(left, left_root, ast.Name(id=left_index)) and _subscript(right, right_root) and _is_index_difference(_slice(right), right_outer, right_inner)) or (_subscript(right, left_root, ast.Name(id=left_index)) and _subscript(left, right_root) and _is_index_difference(_slice(left), right_outer, right_inner))


def _slice(node: ast.AST | None) -> ast.AST | None:
    if not isinstance(node, ast.Subscript):
        return None
    value = node.slice
    if hasattr(ast, "Index") and isinstance(value, ast.Index):
        return value.value
    return value


def _conv_augassign(node: ast.AST, left_root: str, right_root: str, left_index: str, outer: str, inner: str) -> bool:
    return isinstance(node, ast.AugAssign) and isinstance(node.op, ast.Add) and _is_convolution_term(node.value, left_root, left_index, right_root, outer, inner)


def _anchor_record(path: str, node: ast.AST | None, matched: bool, expected: str, reason: str) -> dict[str, Any]:
    return {
        "path": path,
        "matched": bool(matched),
        "expected": expected,
        "reason": reason,
        "line": getattr(node, "lineno", None),
        "ast": _ast_dump(node),
    }


def _anchor_high_mul(tree: ast.AST | None) -> dict[str, Any]:
    node = _method(tree, "IntervalTaylorJet", "__mul__")
    loops = [item for item in ast.walk(node) if isinstance(item, ast.For)] if node else []
    outer = any(_is_range_plus_one(item.iter, "order") for item in loops)
    inner = any(_is_range_plus_one(item.iter, "index") for item in loops)
    convolution = any(_conv_augassign(item, "self", "other", "left", "index", "left") for item in ast.walk(node)) if node else False
    coercion = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) == "self._coerce")
    returned = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) in {"IntervalTaylorJet", "type"})
    ok = _signature(node, ("self", "other")) and outer and inner and convolution and coercion and returned
    return _anchor_record("high.IntervalTaylorJet.__mul__", node, ok, "order/index convolution with self[left]*other[index-left]", "signature, coercion, two loops, convolution, and constructor checked")


def _anchor_high_reciprocal(tree: ast.AST | None) -> dict[str, Any]:
    node = _method(tree, "IntervalTaylorJet", "reciprocal")
    loops = [item for item in ast.walk(node) if isinstance(item, ast.For)] if node else []
    recurrence = any(_conv_augassign(item, "self", "coefficients", "left", "index", "left") for item in ast.walk(node)) if node else False
    nonzero = any(call.args and _subscript(call.args[0], "self") for call in _calls(node, "_require_nonzero"))
    init = _contains(node, lambda item: isinstance(item, ast.BinOp) and isinstance(item.op, ast.Div) and _subscript(item.right, "self"))
    append = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) == "coefficients.append" and isinstance(item.args[0], ast.BinOp) and isinstance(item.args[0].op, ast.Div) and isinstance(item.args[0].left, ast.UnaryOp) and isinstance(item.args[0].left.op, ast.USub) and _subscript(item.args[0].right, "self"))
    loop = any(_is_range_from_one(item.iter, "index") for item in loops)
    returned = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) in {"IntervalTaylorJet", "type"})
    ok = _signature(node, ("self",)) and nonzero and init and recurrence and append and loop and returned
    return _anchor_record("high.IntervalTaylorJet.reciprocal", node, ok, "nonzero self[0] and q_n=-sum(self[k]q[n-k])/self[0]", "denominator guard, recurrence loop, and negative quotient checked")


def _anchor_high_differentiate(tree: ast.AST | None) -> dict[str, Any]:
    node = _method(tree, "IntervalTaylorJet", "differentiate")
    list_comps = [item for item in ast.walk(node) if isinstance(item, ast.ListComp)] if node else []
    matched_comp = False
    for comp in list_comps:
        if len(comp.generators) != 1 or not _name(comp.generators[0].target, "index"):
            continue
        if not isinstance(comp.generators[0].iter, ast.Call) or _dotted(comp.generators[0].iter.func) != "range":
            continue
        expression = comp.elt
        if not isinstance(expression, ast.BinOp) or not isinstance(expression.op, ast.Mult):
            continue
        matched_comp = _is_add_one(expression.left, "index") and _is_add_one(_slice(expression.right), "index")
    returned = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) in {"IntervalTaylorJet", "type"})
    ok = _signature(node, ("self",)) and matched_comp and returned
    return _anchor_record("high.IntervalTaylorJet.differentiate", node, ok, "(index+1)*coefficient[index+1]", "ordinary derivative multiplier checked")


def _anchor_high_source_jet(tree: ast.AST | None) -> dict[str, Any]:
    node = _single_function(tree, "_source_jet")
    owner = "module"
    if node is None:
        owner_class = _class(tree, "IntervalTaylorJet")
        if owner_class is not None:
            candidates = [item for item in owner_class.body if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == "_source_jet"]
            node = candidates[0] if len(candidates) == 1 else None
            owner = "IntervalTaylorJet"
    loops = [item for item in ast.walk(node) if isinstance(item, ast.For)] if node else []
    loop = any(_is_range_plus_one(item.iter, "axial_order") for item in loops)
    source_read = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) == "_source_value" and len(item.args) >= 3 and _name(item.args[1], "radial_order") and _name(item.args[2], "order"))
    factorial = _contains(node, lambda item: _is_factorial_call(item, "order"))
    divided = _contains(node, lambda item: isinstance(item, ast.BinOp) and isinstance(item.op, ast.Div) and _contains(item.right, lambda child: _is_factorial_call(child, "order")))
    returned = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) in {"IntervalTaylorJet", "LaurentTaylorJet", "type"})
    ok = _signature(node, ("context", "grid", "radial_order", "axial_order", "label")) and loop and source_read and factorial and divided and returned
    return _anchor_record("high." + owner + "._source_jet", node, ok, "source derivative rows divided by factorial(order)", "source read, axial loop, factorial division, and return checked")


def _anchor_high_relative(tree: ast.AST | None) -> dict[str, Any]:
    node = _single_function(tree, "_relative_f0_jet")
    loops = [item for item in ast.walk(node) if isinstance(item, ast.For)] if node else []
    outer = any(_is_range_from_one(item.iter, "n") for item in loops)
    inner = any(isinstance(item.iter, ast.Call) and _dotted(item.iter.func) == "range" and len(item.iter.args) == 1 and _name(item.iter.args[0], "n") for item in loops)
    recurrence = any(_contains(item, lambda child: isinstance(child, ast.BinOp) and isinstance(child.op, ast.Mult) and _subscript(child.left, "ell") and _subscript(child.right, "coefficients")) for item in ast.walk(node)) if node else False
    append = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) == "coefficients.append" and item.args and isinstance(item.args[0], ast.BinOp) and isinstance(item.args[0].op, ast.Div) and _name(item.args[0].left, "value") and _is_n_denominator(item.args[0].right))
    ell = _assignments(node, "ell")
    denominator = _assignments(node, "denominator")
    denom_shape_ok = bool(denominator) and _contains(denominator[0].value, lambda item: isinstance(item, ast.BinOp) and isinstance(item.op, ast.Pow) and _integer(item.right, 2) and (_name(item.left, "sigma") or _contains(item.left, lambda child: _name(child, "sigma"))) )
    nonzero = bool(_calls(node, "_require_nonzero"))
    initial_one = _contains(node, lambda item: isinstance(item, ast.List) and item.elts and _integer(item.elts[0], 1)) or _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) == "context.mpf" and item.args and _integer(item.args[0], 1))
    ok = _signature(node, ("context", "z", "delta", "Lambda", "required_j", "order")) and outer and inner and recurrence and append and bool(ell) and bool(denominator) and denom_shape_ok and nonzero and initial_one
    return _anchor_record("high._relative_f0_jet", node, ok, "ell recurrence with denominator H^2+(j/500)^2 and C0=1", "partial AST shapes for recurrence, nonzero guard and sigma exponent checked; full ell and denominator expressions remain assumptions")


def _anchor_low_mul(tree: ast.AST | None) -> dict[str, Any]:
    node = _method(tree, "Dual", "__mul__")
    returns = [item for item in ast.walk(node) if isinstance(item, ast.Call) and _dotted(item.func) == "Dual"] if node else []
    matched = False
    expected_value = ("bin", "*", ("attribute", "self.value"), ("attribute", "other.value"))
    expected_derivative = ("bin", "+", ("bin", "*", ("attribute", "self.derivative"), ("attribute", "other.value")), ("bin", "*", ("attribute", "self.value"), ("attribute", "other.derivative")))
    for call in returns:
        if len(call.args) == 3 and _shape(call.args[1]) == expected_value and _shape(call.args[2]) == expected_derivative:
            matched = True
    coercion = _contains(node, lambda item: isinstance(item, ast.Call) and _dotted(item.func) == "self._coerce")
    ok = _signature(node, ("self", "other")) and matched and coercion
    return _anchor_record("low.Dual.__mul__", node, ok, "Dual(value*value, derivative*value+value*derivative)", "explicit product-rule slots checked")


def _anchor_low_reciprocal(tree: ast.AST | None) -> dict[str, Any]:
    node = _method(tree, "Dual", "reciprocal")
    nonzero = any(call.args and _shape(call.args[0]) == ("attribute", "self.value") for call in _calls(node, "_require_nonzero"))
    inverse = _contains(node, lambda item: isinstance(item, ast.BinOp) and isinstance(item.op, ast.Div) and _shape(item.right) == ("attribute", "self.value"))
    expected = ("unary", "-", ("bin", "/", ("attribute", "self.derivative"), ("bin", "**", ("attribute", "self.value"), ("constant", 2))))
    returned = any(len(call.args) == 3 and _shape(call.args[2]) == expected for call in [item for item in ast.walk(node) if isinstance(item, ast.Call) and _dotted(item.func) == "Dual"]) if node else False
    ok = _signature(node, ("self",)) and nonzero and inverse and returned
    return _anchor_record("low.Dual.reciprocal", node, ok, "Dual(a0^-1,-a1/a0^2)", "denominator guard and squared value in derivative slot checked")


def _anchor_low_fz(tree: ast.AST | None) -> dict[str, Any]:
    node = _single_function(tree, "build")
    assignments = _assignments(node, "fz")
    expected = (
        "bin",
        "+",
        ("bin", "+", ("name", "phi0_2"), ("bin", "*", ("bin", "*", ("constant", 2), ("attribute", "ell.value")), ("name", "phi0_1"))),
        ("bin", "*", ("bin", "+", ("attribute", "ell.derivative"), ("bin", "**", ("attribute", "ell.value"), ("constant", 2))), ("name", "phi0_0")),
    )
    matched = False
    call_node: ast.Call | None = None
    if len(assignments) == 1 and isinstance(assignments[0].value, ast.Call) and _dotted(assignments[0].value.func) == "Dual":
        call_node = assignments[0].value
        matched = len(call_node.args) == 3 and _name(call_node.args[0], "c") and _shape(call_node.args[1]) == ("attribute", "f.derivative") and _shape(call_node.args[2]) == expected
    ok = _signature(node, ("c", "rho", "z", "delta", "Lambda", "j", "Phi_grid", "U_grid", "M_grid")) and matched
    return _anchor_record("low.build.fz", call_node or node, ok, "Dual(c,f.derivative,Phi2+2*ell*Phi1+(ell'+ell^2)*Phi0)", "low fz value and second derivative slot checked")


def _source_receipt(label: str, tree: ast.AST | None, issue: str | None) -> dict[str, Any]:
    return {"label": label, "parsed": tree is not None, "top_level_statement_count": len(tree.body) if isinstance(tree, ast.Module) else None, "error": issue}


def _parse(source: Any, label: str) -> tuple[ast.Module | None, str | None]:
    if not isinstance(source, str):
        return None, "source input must be a string"
    try:
        return ast.parse(source, filename=label, mode="exec"), None
    except SyntaxError as exc:
        return None, f"{exc.msg} at line {exc.lineno}, offset {exc.offset}"


def _all_anchors(high_tree: ast.AST | None, low_tree: ast.AST | None) -> dict[str, dict[str, Any]]:
    return {
        "high_mul": _anchor_high_mul(high_tree),
        "high_reciprocal": _anchor_high_reciprocal(high_tree),
        "high_differentiate": _anchor_high_differentiate(high_tree),
        "high_source_jet": _anchor_high_source_jet(high_tree),
        "high_relative_f0": _anchor_high_relative(high_tree),
        "low_mul": _anchor_low_mul(low_tree),
        "low_reciprocal": _anchor_low_reciprocal(low_tree),
        "low_fz": _anchor_low_fz(low_tree),
    }


def build(high_source: str, low_source: str) -> dict[str, Any]:
    """Parse two source strings and return a conditional exact proof receipt."""

    high_tree, high_issue = _parse(high_source, HIGH_LABEL)
    low_tree, low_issue = _parse(low_source, LOW_LABEL)
    anchors = _all_anchors(high_tree, low_tree)
    anchor_differences = [
        {"kind": "anchor_mismatch", "path": name, "reason": record["reason"]}
        for name, record in anchors.items()
        if not record["matched"]
    ]
    source_differences = []
    if high_issue is not None:
        source_differences.append({"kind": "high_parse", "message": high_issue})
    if low_issue is not None:
        source_differences.append({"kind": "low_parse", "message": low_issue})
    differences = source_differences + anchor_differences
    proof = _exact_algebra_report()
    checks = {
        "high_parsed": high_tree is not None,
        "low_parsed": low_tree is not None,
        "all_ast_anchors": not anchor_differences,
        "exact_scalar_theorem": proof["exact"],
    }
    passed_count = sum(1 for value in checks.values() if value)
    passed = all(checks.values()) and not differences
    if passed:
        status = "proved"
    elif high_tree is not None or low_tree is not None:
        status = "partial"
    else:
        status = "rejected"
    return {
        "schema": SCHEMA,
        "status": status,
        "passed": passed,
        "ok": passed,
        "conditional": True,
        "claim": "exact first-order Taylor-ring projection to the low Dual formulas",
        "source": {
            "high": _source_receipt(HIGH_LABEL, high_tree, high_issue),
            "low": _source_receipt(LOW_LABEL, low_tree, low_issue),
        },
        "anchors": anchors,
        "checks": {"passed": passed_count, "total": len(checks), "items": checks},
        "theorem": proof,
        "exact_computational_theorem": proof,
        "assumptions": list(ASSUMPTIONS),
        "unproved": list(UNPROVED),
        "scope_notes": list(SCOPE_NOTES),
        "hard_to_certify": list(HARD_TO_CERTIFY),
        "physical_field": False,
        "sum_claim": False,
        "physical_field_materialized": False,
        "sum_over_fields_claimed": False,
        "differences": differences,
    }


__all__ = [
    "ASSUMPTIONS",
    "HIGH_LABEL",
    "LOW_LABEL",
    "SCHEMA",
    "HARD_TO_CERTIFY",
    "SCOPE_NOTES",
    "UNPROVED",
    "SourceIdentityError",
    "build",
]
