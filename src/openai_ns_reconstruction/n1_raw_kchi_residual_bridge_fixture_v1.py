"""Exact Fraction witnesses for :mod:`n1_raw_kchi_residual_bridge_v1`.

The fixture is intentionally self-contained.  It uses no historical module,
no floating point, no ``mpmath``, and no interval conversion.  Its private
2-by-2 and 3-by-3 matrix/vector classes perform only exact ``Fraction``
arithmetic.  The matrices are non-commuting, and the ``A_Z`` matrices are
nonzero so that the ordinary-Z1 transport cannot pass by accident.
The last-row/first-column entries are the finite-dimensional E51 analogues
(E21 and E31 in the two synthetic dimensions); the bridge contract retains
the canonical ``A=-2*chi*E51`` and ``A_Z=-2*chi_Z*E51`` symbols.

``run_fixture()`` is a caller-run check rather than an import-time test.  For
each size it hand-forms ``QW`` and ``Q(W_Z)`` (including ``B1*W_ZZ``),
``d/dZ(QW)`` with nonzero ``C0_Z*W`` and ``B1_Z*W_Z``, then
``d_raw``, ``d_raw_Z``, ``d_chi``, and ``d_chi_Z``.  It then compares the
bridge output with those direct residual functionals and checks exact zero
for the two bridge identities.  Negative cases cover a wrong A sign, owner
mismatch, missing certificate, unverified resolvent, nonzero/unknown ``G_Z``,
missing Z1, incomplete Q derivative certification, and mismatched units.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Iterable

from .n1_raw_kchi_residual_bridge_v1 import (
    A_DEFINITION,
    AZ_DEFINITION,
    DQW_D_Z_DEFINITION,
    Q_ON_W_Z_DEFINITION,
    Q_DEFINITION,
    ResidualBridgeError,
    bridge,
)


SCHEMA = "n1_raw_kchi_residual_bridge_fixture_v1"


def _f(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    return Fraction(value)


class _Vector:
    """Minimal exact vector field with the caller-required ``+`` operation."""

    __slots__ = ("values",)

    def __init__(self, values: Iterable[Any]) -> None:
        self.values = tuple(_f(value) for value in values)

    def __add__(self, other: Any) -> "_Vector":
        if not isinstance(other, _Vector) or len(self.values) != len(other.values):
            return NotImplemented
        return _Vector(left + right for left, right in zip(self.values, other.values))

    def __sub__(self, other: Any) -> "_Vector":
        if not isinstance(other, _Vector) or len(self.values) != len(other.values):
            return NotImplemented
        return _Vector(left - right for left, right in zip(self.values, other.values))

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, _Vector) and self.values == other.values

    def __repr__(self) -> str:
        return "Vector(" + ", ".join(str(value) for value in self.values) + ")"


class _Matrix:
    """Small exact square matrix used only by this synthetic fixture."""

    __slots__ = ("rows", "size")

    def __init__(self, rows: Iterable[Iterable[Any]]) -> None:
        normalized = tuple(tuple(_f(value) for value in row) for row in rows)
        if not normalized or any(len(row) != len(normalized) for row in normalized):
            raise ValueError("fixture matrices must be nonempty and square")
        self.rows = normalized
        self.size = len(normalized)

    def __matmul__(self, other: Any) -> Any:
        if isinstance(other, _Vector):
            if len(other.values) != self.size:
                raise ValueError("matrix/vector dimensions do not match")
            return _Vector(
                sum((self.rows[row][column] * other.values[column]
                     for column in range(self.size)), Fraction(0))
                for row in range(self.size)
            )
        if isinstance(other, _Matrix):
            if other.size != self.size:
                raise ValueError("matrix dimensions do not match")
            return _Matrix(
                tuple(
                    sum((self.rows[row][inner] * other.rows[inner][column]
                         for inner in range(self.size)), Fraction(0))
                    for column in range(self.size)
                )
                for row in range(self.size)
            )
        return NotImplemented

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, _Matrix) and self.rows == other.rows

    def __repr__(self) -> str:
        return "Matrix(" + repr(self.rows) + ")"


class _LinearOperator:
    """Callable operator carrying an explicit owner and optional A formula."""

    def __init__(self, matrix: _Matrix, owner: Any, *, definition: str | None = None) -> None:
        self.matrix = matrix
        self.owner = owner
        if definition is not None:
            self.definition = definition

    def __call__(self, field: _Vector) -> _Vector:
        result = self.matrix @ field
        if not isinstance(result, _Vector):
            raise TypeError("fixture operator expects a vector field")
        return result


class _Backend:
    """Backend path witness; it exposes operators plus field addition."""

    def __init__(self, owner: Any, operators: dict[str, _LinearOperator]) -> None:
        self.owner = owner
        for name, operator in operators.items():
            setattr(self, name, operator)

    @staticmethod
    def add(left: _Vector, right: _Vector) -> _Vector:
        return left + right


def _identity(size: int) -> _Matrix:
    return _Matrix(
        tuple(
            tuple(Fraction(1) if row == column else Fraction(0) for column in range(size))
            for row in range(size)
        )
    )


def _certificate(owner: Any, *, G_Z: Any = 0, resolvent: bool = True) -> dict[str, Any]:
    return {
        "owner": owner,
        "resolvent_relation_verified": resolvent,
        "G_Z": G_Z,
        "z1_provided": True,
        "A_definition": A_DEFINITION,
        "Az_definition": AZ_DEFINITION,
        "q_definition": Q_DEFINITION,
        "q_on_w_z_definition": Q_ON_W_Z_DEFINITION,
        "d_qw_d_z_definition": DQW_D_Z_DEFINITION,
        "q_w_z_includes_B1_W_ZZ": True,
        "d_qw_d_z_includes_C0_Z_W": True,
        "d_qw_d_z_includes_B1_Z_W_Z": True,
        "d_raw_Z_includes_full_dQW": True,
        "units": "synthetic-field-units",
        "coordinates": "fixed-x,Z",
        "initial_data": "same-W-and-g",
        "same_units": True,
        "same_coordinates": True,
        "same_initial_data": True,
    }


def _direct_case(
    *,
    G: _Matrix,
    K: _Matrix,
    A: _Matrix,
    Az: _Matrix,
    C0: _Matrix,
    C0_Z: _Matrix,
    B1: _Matrix,
    B1_Z: _Matrix,
    W: _Vector,
    W_Z: _Vector,
    W_ZZ: _Vector,
    g: _Vector,
    g_Z: _Vector,
) -> dict[str, Any]:
    """Hand-form every residual row for one exact finite-dimensional witness."""

    QW = (C0 @ W) + (B1 @ W_Z)
    # Q(W_Z) is not the same object as d/dZ(QW) when C0 or B1 depends on Z.
    Q_on_W_Z = (C0 @ W_Z) + (B1 @ W_ZZ)
    dQW_dZ = (C0_Z @ W) + Q_on_W_Z + (B1_Z @ W_Z)
    d_raw = W - (G @ (g + QW))
    d_raw_Z = W_Z - (G @ (g_Z + dQW_dZ))
    chi_input = g + QW - (A @ W)
    chi_input_Z = g_Z + dQW_dZ - (Az @ W) - (A @ W_Z)
    d_chi = W - (K @ chi_input)
    # K_Z = K*A_Z*K follows from K=(I-G*A)^(-1)G with fixed G.  Keeping
    # this term is the direct residual derivative that the bridge identity
    # transports through G*A_Z*d_chi.
    K_Z = K @ Az @ K
    d_chi_Z = W_Z - (K_Z @ chi_input) - (K @ chi_input_Z)
    return {
        "QW": QW,
        "Q_on_W_Z": Q_on_W_Z,
        "dQW_dZ": dQW_dZ,
        "C0_Z_W": C0_Z @ W,
        "B1_Z_W_Z": B1_Z @ W_Z,
        "d_raw": d_raw,
        "d_raw_Z": d_raw_Z,
        "d_chi": d_chi,
        "d_chi_Z": d_chi_Z,
    }


def _case_2x2(owner: Any) -> dict[str, Any]:
    """A noncommuting 2-by-2 exact witness with chi_Z = 1/4."""

    G = _Matrix(((1, 2), (0, 1)))
    # chi=1/3, so A=-2*chi*E21.  The hard-coded K is
    # (I-GA)^(-1)G, not G(I-AG)^(-1); the order is material.
    A = _Matrix(((0, 0), (Fraction(-2, 3), 0)))
    Az = _Matrix(((0, 0), (Fraction(-1, 2), 0)))  # chi_Z=1/4
    K = _Matrix(((Fraction(3, 7), Fraction(6, 7)),
                 (Fraction(-2, 7), Fraction(3, 7))))
    C0 = _Matrix(((1, 1), (0, 2)))
    C0_Z = _Matrix(((1, 0), (0, 0)))
    B1 = _Matrix(((0, 1), (1, 0)))
    B1_Z = _Matrix(((0, 0), (1, 0)))
    W = _Vector((2, -1))
    W_Z = _Vector((1, 3))
    W_ZZ = _Vector((2, -1))
    g = _Vector((Fraction(1, 3), Fraction(2, 5)))
    g_Z = _Vector((Fraction(-1, 2), 1))
    direct = _direct_case(
        G=G, K=K, A=A, Az=Az, C0=C0, C0_Z=C0_Z, B1=B1, B1_Z=B1_Z,
        W=W, W_Z=W_Z, W_ZZ=W_ZZ, g=g, g_Z=g_Z,
    )
    operators = {
        "G": _LinearOperator(G, owner),
        "K": _LinearOperator(K, owner),
        "A": _LinearOperator(A, owner, definition=A_DEFINITION),
        "Az": _LinearOperator(Az, owner, definition=AZ_DEFINITION),
    }
    return {
        "dimension": 2,
        "G": G,
        "K": K,
        "A": A,
        "Az": Az,
        "C0": C0,
        "C0_Z": C0_Z,
        "B1": B1,
        "B1_Z": B1_Z,
        "W": W,
        "W_Z": W_Z,
        "W_ZZ": W_ZZ,
        "g": g,
        "g_Z": g_Z,
        "direct": direct,
        "operators": operators,
        "certificate": _certificate(owner),
        "chi_Z_nonzero": True,
    }


def _case_3x3(owner: Any) -> dict[str, Any]:
    """A noncommuting 3-by-3 exact witness with chi_Z = 1/7."""

    G = _Matrix(((1, 2, 0), (0, 1, 1), (0, 0, 1)))
    # chi=1/5, so A=-2*chi*E31.  Again K is hand-expanded from the
    # left resolvent (I-GA)^(-1)G.
    A = _Matrix(((0, 0, 0), (0, 0, 0), (Fraction(-2, 5), 0, 0)))
    Az = _Matrix(((0, 0, 0), (0, 0, 0), (Fraction(-2, 7), 0, 0)))  # chi_Z=1/7
    K = _Matrix(((1, 2, 0),
                 (Fraction(-2, 5), Fraction(1, 5), 1),
                 (Fraction(-2, 5), Fraction(-4, 5), 1)))
    C0 = _Matrix(((1, 1, 0), (0, 2, 1), (1, 0, 1)))
    C0_Z = _Matrix(((0, 1, 0), (0, 0, 0), (0, 0, 0)))
    B1 = _Matrix(((0, 1, 1), (1, 0, 0), (0, 1, 0)))
    B1_Z = _Matrix(((0, 0, 0), (1, 0, 0), (0, 0, 0)))
    W = _Vector((1, -2, 3))
    W_Z = _Vector((2, 1, -1))
    W_ZZ = _Vector((-1, 2, 1))
    g = _Vector((Fraction(1, 2), Fraction(-1, 3), Fraction(2, 5)))
    g_Z = _Vector((Fraction(1, 4), Fraction(2, 3), Fraction(-1, 2)))
    direct = _direct_case(
        G=G, K=K, A=A, Az=Az, C0=C0, C0_Z=C0_Z, B1=B1, B1_Z=B1_Z,
        W=W, W_Z=W_Z, W_ZZ=W_ZZ, g=g, g_Z=g_Z,
    )
    operators = {
        "G": _LinearOperator(G, owner),
        "K": _LinearOperator(K, owner),
        "A": _LinearOperator(A, owner, definition=A_DEFINITION),
        "Az": _LinearOperator(Az, owner, definition=AZ_DEFINITION),
    }
    return {
        "dimension": 3,
        "G": G,
        "K": K,
        "A": A,
        "Az": Az,
        "C0": C0,
        "C0_Z": C0_Z,
        "B1": B1,
        "B1_Z": B1_Z,
        "W": W,
        "W_Z": W_Z,
        "W_ZZ": W_ZZ,
        "g": g,
        "g_Z": g_Z,
        "direct": direct,
        "operators": operators,
        "certificate": _certificate(owner),
        "chi_Z_nonzero": True,
    }


def _rejected(call: Any) -> bool:
    try:
        call()
    except ResidualBridgeError:
        return True
    return False


def _check_case(case: dict[str, Any], *, use_backend: bool) -> dict[str, bool]:
    owner = case["certificate"]["owner"]
    operators = case["operators"]
    backend = _Backend(owner, operators) if use_backend else None
    explicit = {} if use_backend else operators
    direct = case["direct"]
    result = bridge(
        direct["d_raw"],
        direct["d_raw_Z"],
        owner=owner,
        certificate=case["certificate"],
        backend=backend,
        **explicit,
    )
    omitted_B1_WZZ = case["C0"] @ case["W_Z"]
    zero = _Vector((Fraction(0),) * case["dimension"])
    checks = {
        "direct_d_chi_exact": result["d_chi"] == direct["d_chi"],
        "direct_d_chi_Z_exact": result["d_chi_Z"] == direct["d_chi_Z"],
        "d_chi_bridge_exact_zero": (
            result["d_chi"] - (direct["d_raw"] + (case["K"] @ (case["A"] @ direct["d_raw"])))
            == _Vector((Fraction(0),) * case["dimension"])
        ),
        "d_chi_Z_bridge_exact_zero": (
            result["d_chi_Z"]
            - (
                direct["d_raw_Z"]
                + (case["G"] @ (case["Az"] @ result["d_chi"]))
                + (
                    case["K"]
                    @ (case["A"] @ (
                        direct["d_raw_Z"]
                        + (case["G"] @ (case["Az"] @ result["d_chi"]))
                    ))
                )
            )
            == _Vector((Fraction(0),) * case["dimension"])
        ),
        "noncommuting_GA": case["G"] @ case["A"] != case["A"] @ case["G"],
        "chi_Z_nonzero": case["chi_Z_nonzero"],
        "Q_on_W_Z_has_B1_W_ZZ": direct["Q_on_W_Z"] != omitted_B1_WZZ,
        "dQW_dZ_has_C0_Z_W": direct["C0_Z_W"] != zero,
        "dQW_dZ_has_B1_Z_W_Z": direct["B1_Z_W_Z"] != zero,
        "Q_on_W_Z_distinct_from_full_dQW": direct["Q_on_W_Z"] != direct["dQW_dZ"],
        "d_raw_Z_uses_full_dQW": (
            direct["d_raw_Z"]
            == case["W_Z"] - (case["G"] @ (case["g_Z"] + direct["dQW_dZ"]))
        ),
        "functional_output_not_bound": (
            result["source_functional"]["is_value"] is False
            and result["source_functional"]["is_error_bound"] is False
            and result["signed_residual_contract"]["zero_assumed"] is False
        ),
        "certificate_is_premise": (
            result["certificate"]["certificate_is_caller_premise"] is True
            and result["certificate"]["actual_source_verified"] is False
        ),
    }
    return checks


def run_fixture() -> dict[str, Any]:
    """Run the exact synthetic checks and return a serializable check ledger."""

    owner_2 = object()
    owner_3 = object()
    case_2 = _case_2x2(owner_2)
    case_3 = _case_3x3(owner_3)
    checks: dict[str, bool] = {}
    checks.update({f"2x2_{key}": value for key, value in _check_case(case_2, use_backend=False).items()})
    checks.update({f"3x3_{key}": value for key, value in _check_case(case_3, use_backend=True).items()})

    # Keep the negative cases explicit and owner-identity based.  A wrong sign
    # is carried by the operator itself, so a matching-looking label cannot
    # make it admissible.
    wrong_A = _LinearOperator(
        _Matrix(((0, 0), (Fraction(2, 3), 0))),
        owner_2,
        definition="A=+2*chi*E21",
    )
    checks["wrong_sign_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=case_2["certificate"],
            G=case_2["operators"]["G"],
            K=case_2["operators"]["K"],
            A=wrong_A,
            Az=case_2["operators"]["Az"],
        )
    )
    checks["owner_mismatch_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=object(),
            certificate=case_2["certificate"],
            **case_2["operators"],
        )
    )
    checks["missing_certificate_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            **case_2["operators"],
        )
    )
    no_resolvent = dict(case_2["certificate"], resolvent_relation_verified=False)
    checks["unverified_resolvent_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=no_resolvent,
            **case_2["operators"],
        )
    )
    bad_gz = dict(case_2["certificate"], G_Z=1)
    checks["nonzero_G_Z_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=bad_gz,
            **case_2["operators"],
        )
    )
    float_gz = dict(case_2["certificate"], G_Z=0.0)
    checks["float_zero_G_Z_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=float_gz,
            **case_2["operators"],
        )
    )
    unknown_gz = dict(case_2["certificate"])
    unknown_gz.pop("G_Z")
    checks["unknown_G_Z_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=unknown_gz,
            **case_2["operators"],
        )
    )
    checks["missing_Z1_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            None,
            owner=owner_2,
            certificate=case_2["certificate"],
            **case_2["operators"],
        )
    )
    mismatched_units = dict(case_2["certificate"], raw_units="u", chi_units="v")
    checks["different_units_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=mismatched_units,
            **case_2["operators"],
        )
    )
    missing_context = dict(case_2["certificate"])
    for key in (
        "units",
        "coordinates",
        "initial_data",
        "same_units",
        "same_coordinates",
        "same_initial_data",
    ):
        missing_context.pop(key, None)
    checks["missing_context_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=missing_context,
            **case_2["operators"],
        )
    )
    missing_q = dict(case_2["certificate"])
    for key in (
        "q_definition",
        "q_on_w_z_definition",
        "d_qw_d_z_definition",
        "q_w_z_includes_B1_W_ZZ",
        "d_qw_d_z_includes_C0_Z_W",
        "d_qw_d_z_includes_B1_Z_W_Z",
        "d_raw_Z_includes_full_dQW",
    ):
        missing_q.pop(key, None)
    checks["missing_Q_contract_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=missing_q,
            **case_2["operators"],
        )
    )
    incomplete_q = dict(case_2["certificate"], d_qw_d_z_includes_B1_Z_W_Z=False)
    checks["incomplete_Q_derivative_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=incomplete_q,
            **case_2["operators"],
        )
    )
    bad_az_definition = dict(case_2["certificate"], Az_definition=A_DEFINITION)
    checks["wrong_A_Z_definition_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=bad_az_definition,
            **case_2["operators"],
        )
    )
    checks["ownerless_operator_rejected"] = _rejected(
        lambda: bridge(
            case_2["direct"]["d_raw"],
            case_2["direct"]["d_raw_Z"],
            owner=owner_2,
            certificate=case_2["certificate"],
            G=lambda field: field,
            K=lambda field: field,
            A=lambda field: field,
            Az=lambda field: field,
        )
    )
    checks["identity_matrices_exact"] = (
        case_2["G"] @ case_2["A"] != case_2["A"] @ case_2["G"]
        and case_3["G"] @ case_3["A"] != case_3["A"] @ case_3["G"]
    )
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "dimensions": [2, 3],
        "arithmetic": "exact Fraction only",
        "source_modules_called": ["n1_raw_kchi_residual_bridge_v1"],
        "historical_modules_called": False,
        "mpmath_used": False,
        "signed_residual_zero_assumed": False,
        "norm_cap_as_field": False,
    }


def run() -> dict[str, Any]:
    """Short alias for :func:`run_fixture`."""

    return run_fixture()


__all__ = ["SCHEMA", "run_fixture", "run"]

