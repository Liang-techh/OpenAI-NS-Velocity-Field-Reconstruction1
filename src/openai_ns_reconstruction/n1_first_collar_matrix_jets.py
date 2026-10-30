"""Second-order ordinary-``Z`` jets for the raw n=1 first-collar matrix.

This is a small, dependency-free numerical leaf.  It transcribes the two
symbolic dictionaries returned by ``raw_matrix_expressions()`` in
``experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py``.
The source is cited as provenance only: this module does not import it, read
it at run time, or consume an old result/JSON packet.

The public evaluator accepts scalar ``Z``-constant parameters ``rho``, ``x``,
``delta`` and ``epsilon``; a scalar ``Z`` value (whose ordinary derivative is
set to one); and explicit ordinary, *unfactorialised* three-term jets for
``F``, ``Fr``, ``FZ``, ``U``, ``Ur``, ``UZ`` and ``Q``.  A jet is ordered as
``(value, d/dZ, d^2/dZ^2)``.  The result is a mapping with exactly ``C0`` and
``B1`` keys.  Each matrix entry is itself a three-term ordinary jet tuple,
and matrix keys are the signed row/column labels from the source (for
example ``"5,2"``).

The implementation deliberately never extracts interval endpoints or
midpoints.  Values such as :class:`mpmath.iv.ivmpf` therefore remain interval
objects through every operation; only the four explicitly constant
parameters have zero ``Z`` derivatives.  Division requires the denominator's
value component to be invertible in the supplied numeric type, as it does in
the source expressions (notably ``L = 1 - delta*Z**2``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


COMPONENT = "n1_first_collar_matrix_jets"
ORDERS = (0, 1, 2)
MATRIX_KEYS = {
    "C0": (
        "4,1",
        "5,1",
        "5,2",
        "5,3",
        "5,5",
        "6,1",
        "6,2",
        "6,3",
        "6,4",
        "6,6",
    ),
    "B1": ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4"),
}


@dataclass(frozen=True, slots=True)
class Jet2:
    """An ordinary (not factorial-divided) second-order ``Z`` jet.

    ``value`` is the zeroth derivative, ``first`` is the first derivative,
    and ``second`` is the second derivative.  Components intentionally have
    no prescribed scalar type: Python integers, ``decimal`` values,
    ``mpmath`` values, and ``mpmath.iv`` intervals all work as long as their
    arithmetic supports the requested operations.
    """

    value: Any
    first: Any
    second: Any

    @property
    def d0(self) -> Any:
        """Alias for the zeroth ordinary derivative."""

        return self.value

    @property
    def d1(self) -> Any:
        """Alias for the first ordinary derivative."""

        return self.first

    @property
    def d2(self) -> Any:
        """Alias for the second ordinary derivative."""

        return self.second

    def __iter__(self):
        return iter((self.value, self.first, self.second))

    def __len__(self) -> int:
        return 3

    def __getitem__(self, index: int) -> Any:
        return (self.value, self.first, self.second)[index]

    def __add__(self, other: Any) -> "Jet2":
        return jet_add(self, other)

    def __radd__(self, other: Any) -> "Jet2":
        return jet_add(other, self)

    def __sub__(self, other: Any) -> "Jet2":
        return jet_sub(self, other)

    def __rsub__(self, other: Any) -> "Jet2":
        return jet_sub(other, self)

    def __mul__(self, other: Any) -> "Jet2":
        return jet_mul(self, other)

    def __rmul__(self, other: Any) -> "Jet2":
        return jet_mul(other, self)

    def __truediv__(self, other: Any) -> "Jet2":
        return jet_div(self, other)

    def __rtruediv__(self, other: Any) -> "Jet2":
        return jet_div(other, self)

    def __neg__(self) -> "Jet2":
        return Jet2(-self.value, -self.first, -self.second)


def jet_constant(value: Any) -> Jet2:
    """Lift a scalar to a ``Z``-constant jet."""

    return Jet2(value, 0, 0)


def jet_variable(value: Any) -> Jet2:
    """Lift a scalar ``Z`` value to ``(value, 1, 0)``."""

    return Jet2(value, 1, 0)


def _as_jet(value: Any) -> Jet2:
    if isinstance(value, Jet2):
        return value
    return jet_constant(value)


def jet_add(*terms: Any) -> Jet2:
    """Add zero or more scalar/jet operands componentwise."""

    if not terms:
        return jet_constant(0)
    result = _as_jet(terms[0])
    for term in terms[1:]:
        right = _as_jet(term)
        result = Jet2(
            result.value + right.value,
            result.first + right.first,
            result.second + right.second,
        )
    return result


def jet_sub(left: Any, right: Any) -> Jet2:
    """Subtract two scalar/jet operands componentwise."""

    a, b = _as_jet(left), _as_jet(right)
    return Jet2(a.value - b.value, a.first - b.first, a.second - b.second)


def jet_mul(left: Any, right: Any) -> Jet2:
    """Multiply two scalar/jet operands using the ordinary product rule."""

    a, b = _as_jet(left), _as_jet(right)
    return Jet2(
        a.value * b.value,
        a.first * b.value + a.value * b.first,
        a.second * b.value + 2 * a.first * b.first + a.value * b.second,
    )


def jet_inv(value: Any) -> Jet2:
    """Invert a scalar/jet using ordinary first and second derivatives."""

    a = _as_jet(value)
    a0, a1, a2 = a.value, a.first, a.second
    inv0 = 1 / a0
    return Jet2(
        inv0,
        -a1 / (a0 * a0),
        2 * a1 * a1 / (a0 * a0 * a0) - a2 / (a0 * a0),
    )


def jet_div(left: Any, right: Any) -> Jet2:
    """Divide two scalar/jet operands using :func:`jet_inv`."""

    return jet_mul(left, jet_inv(right))


def _parameter(value: Any, name: str) -> Jet2:
    """Require and lift one of the four scalar ``Z``-constant parameters."""

    if isinstance(value, Jet2):
        raise TypeError(f"{name} must be a scalar Z-constant, not a jet")
    return jet_constant(value)


def _source_jet(value: Any, name: str) -> Jet2:
    """Coerce one required explicit three-component source jet."""

    if isinstance(value, Jet2):
        return value
    if isinstance(value, (tuple, list)) and len(value) == 3:
        return Jet2(value[0], value[1], value[2])
    raise TypeError(f"{name} must be an ordinary Z-jet with three components")


def _z_jet(value: Any) -> Jet2:
    """Coerce the independent ``Z`` value, whose default derivative is one.

    The documented form is a scalar.  Accepting an explicit :class:`Jet2` or
    a length-two/length-three sequence is useful when a caller has already
    materialised the independent variable; in those forms the caller is
    responsible for supplying first derivative one and second derivative
    zero.
    """

    if isinstance(value, Jet2):
        return value
    if isinstance(value, (tuple, list)):
        if len(value) == 2:
            return Jet2(value[0], value[1], 0)
        if len(value) == 3:
            return Jet2(value[0], value[1], value[2])
        raise TypeError("Z sequence must contain two or three components")
    return jet_variable(value)


def _entry_tuple(value: Jet2) -> tuple[Any, Any, Any]:
    """Expose a jet as the public ordinary derivative tuple."""

    return (value.value, value.first, value.second)


def evaluate_raw_matrices(
    rho: Any,
    x: Any,
    delta: Any,
    epsilon: Any,
    Z: Any,
    F: Any,
    Fr: Any,
    FZ: Any,
    U: Any,
    Ur: Any,
    UZ: Any,
    Q: Any,
) -> dict[str, dict[str, tuple[Any, Any, Any]]]:
    """Evaluate the raw first-collar ``C0`` and ``B1`` matrix jets.

    Parameters ``rho``, ``x``, ``delta`` and ``epsilon`` are scalar constants
    with respect to ``Z``.  ``Z`` is normally a scalar value and is lifted to
    ``(Z, 1, 0)``.  ``F``, ``Fr``, ``FZ``, ``U``, ``Ur``, ``UZ`` and ``Q`` are
    each explicit three-component ordinary jets, with no factorial division.

    The returned mapping has the shape ``{"C0": {"row,col": (j0,j1,j2)},
    "B1": {...}}``.  ``C0`` is exactly source ``B0`` after removing the
    fixed kinematic entries ``1,5``, ``2,6`` and ``3,6``.
    """

    r = _parameter(rho, "rho")
    xx = _parameter(x, "x")
    de = _parameter(delta, "delta")
    eps = _parameter(epsilon, "epsilon")
    z = _z_jet(Z)
    f = _source_jet(F, "F")
    fr = _source_jet(Fr, "Fr")
    fz = _source_jet(FZ, "FZ")
    u = _source_jet(U, "U")
    ur = _source_jet(Ur, "Ur")
    uz = _source_jet(UZ, "UZ")
    q = _source_jet(Q, "Q")

    one = jet_constant(1)
    two = jet_constant(2)

    # These names and parenthesisations follow raw_matrix_expressions() and
    # intentionally remain visible so that each source term can be audited.
    d = one - z * z
    L = one - de * z * z
    H = jet_div(z * (one - de), two) + d * u
    B = q + (one - two * z * u) / L
    E = z * u - jet_div(one, two)
    Cf = f + r * fr
    Zf = ((-2 - de) * z * f + d * fz - 2 * z * r * fr) / L
    Zw = ((-1 - de) * z * u + d * uz - 2 * z * r * ur) / L

    B0: dict[str, Jet2] = {
        "4,1": eps * 4 * xx * f,
        "5,1": eps * 2 * (q + E * (-2 + de) / L),
        "5,2": eps * (2 * Zf - 2 * (-1 + de) * z * Cf / L),
        "5,3": -eps * 2 * (1 + de) * z * Cf / L,
        "5,5": eps * xx * B,
        "6,1": -eps * eps * 8 * z * r * f / L,
        "6,2": eps * 2 * (E * (-1 + de) / L + Zw - ur * z * (-1 + de) * r / L),
        "6,3": -eps * 2 * (1 + de) * z * r * ur / L,
        "6,4": -eps * 4 * z / L,
        "6,6": eps * xx * B,
    }
    B1: dict[str, Jet2] = {
        "5,1": eps * 2 * H / L,
        "5,2": -eps * 2 * d * Cf / L,
        "5,3": -eps * 2 * d * Cf / L,
        "6,2": eps * 2 * (H - ur * d * r) / L,
        "6,3": -eps * 2 * d * r * ur / L,
        "6,4": eps * 2 * d / L,
    }

    # Keep the output strictly to the requested C0/B1 signed entries and
    # return ordinary derivative triples rather than interval endpoints.
    return {
        "C0": {key: _entry_tuple(B0[key]) for key in MATRIX_KEYS["C0"]},
        "B1": {key: _entry_tuple(B1[key]) for key in MATRIX_KEYS["B1"]},
    }


# Descriptive aliases keep the leaf convenient for callers that name the
# operation after the source function or after the first collar itself.
raw_collar_matrix_jets = evaluate_raw_matrices
first_collar_matrix_jets = evaluate_raw_matrices
raw_matrix_jets = evaluate_raw_matrices


__all__ = [
    "COMPONENT",
    "ORDERS",
    "MATRIX_KEYS",
    "Jet2",
    "jet_constant",
    "jet_variable",
    "jet_add",
    "jet_sub",
    "jet_mul",
    "jet_inv",
    "jet_div",
    "evaluate_raw_matrices",
    "raw_collar_matrix_jets",
    "first_collar_matrix_jets",
    "raw_matrix_jets",
]
