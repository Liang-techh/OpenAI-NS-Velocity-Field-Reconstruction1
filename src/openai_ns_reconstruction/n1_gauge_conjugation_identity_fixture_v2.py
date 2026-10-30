"""Independent exact fixture for :mod:`n1_gauge_conjugation_identity_v1`.

The fixture uses only small ``Fraction`` matrices and independently writes
the product formulas.  It exercises a nonconstant diagonal ``S``, the
explicit six-component natural extra gauge (including an independent
``Pstar_squared`` fourth entry), a state and its ordinary-Z derivative, and a
nonzero ``chi``/``chi_Z`` pair.  It is a synthetic algebra check only: no
source owner, receipt, JSON artifact, historical builder, or physical field
is consulted.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable, Mapping


SIZE = 6
E51 = (4, 0)
REAL_EXTRA_FORMULA = (
    "1", "Lambda^-1", "Lambda^-1", "Pstar_squared*Lambda^-1",
    "1", "Lambda^-1",
)


def _matmul(left: list[list[Any]], right: list[list[Any]]) -> list[list[Any]]:
    return [
        [sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction(0))
         for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def _add(left: list[list[Any]], right: list[list[Any]], sign: int = 1) -> list[list[Any]]:
    return [[a + sign * b for a, b in zip(row_a, row_b)]
            for row_a, row_b in zip(left, right)]


def _diag(entries: list[Any]) -> list[list[Any]]:
    return [[entries[i] if i == j else Fraction(0) for j in range(len(entries))]
            for i in range(len(entries))]


def _inv(entries: list[Any]) -> list[list[Any]]:
    return _diag([1 / value for value in entries])


def _vec(matrix: list[list[Any]], vector: list[Any]) -> list[Any]:
    return [sum((matrix[i][j] * vector[j] for j in range(len(vector))), Fraction(0))
            for i in range(len(matrix))]


def _same(left: Any, right: Any) -> bool:
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            return False
        return all(_same(a, b) for a, b in zip(left, right))
    return left == right


def synthetic_inputs() -> dict[str, Any]:
    """Return independent exact inputs with nonzero off-diagonal raw rows."""

    q = Fraction
    a0 = [[q(0) for _ in range(SIZE)] for _ in range(SIZE)]
    a1 = [[q(0) for _ in range(SIZE)] for _ in range(SIZE)]
    for i in range(SIZE):
        a0[i][i] = q(i + 2, 7)
        a1[i][i] = q(2 * i + 1, 11)
    for i, j, value in ((0, 1, q(3, 5)), (1, 4, q(-2, 7)),
                        (2, 0, q(5, 13)), (4, 3, q(-4, 9)),
                        (5, 2, q(7, 8))):
        a0[i][j] = value
    for i, j, value in ((0, 5, q(-1, 3)), (1, 2, q(4, 7)),
                        (3, 0, q(2, 5)), (4, 1, q(-3, 10)),
                        (5, 4, q(6, 11))):
        a1[i][j] = value

    s = [q(2), q(3), q(5), q(7), q(2), q(11)]
    s_z = [q(1, 3), q(-2, 5), q(3, 7), q(4, 9), q(-1, 6), q(5, 8)]
    akin = [[q(0) for _ in range(SIZE)] for _ in range(SIZE)]
    akin[0][4] = q(1)
    akin[1][5] = q(2)
    akin[2][5] = q(-3)
    chi = q(3, 4)
    chi_z = q(-2, 9)
    achi = [[q(0) for _ in range(SIZE)] for _ in range(SIZE)]
    achi[E51[0]][E51[1]] = -2 * chi
    return {
        "A0": a0,
        "A1": a1,
        "S": _diag(s),
        "S_Z": _diag(s_z),
        "Akin": akin,
        "Achi": achi,
        # The real natural-component extra gauge is supplied explicitly.
        # Pstar_squared is independent of epsilon and is intentionally kept
        # separate in the synthetic data rather than folded into a scalar G.
        "epsilon": q(2, 5),
        "Pstar_squared": q(7, 3),
        "T": _diag([q(1), q(2, 5), q(2, 5), q(14, 15), q(1), q(2, 5)]),
        "T_Z": _diag([q(0) for _ in range(SIZE)]),
        "extra_descriptor": {
            "kind": "natural_component",
            "formula": REAL_EXTRA_FORMULA,
            "Pstar_squared_source": "synthetic-independent-parameter",
            "Pstar_squared_value": q(7, 3),
        },
        "chi": chi,
        "chi_Z": chi_z,
        "Y": [q(1, 2), q(-2, 3), q(3, 5), q(4, 7), q(-5, 11), q(7, 13)],
        "Y_Z": [q(-1, 4), q(2, 5), q(-3, 7), q(5, 8), q(6, 11), q(-7, 10)],
    }


def _reference(data: Mapping[str, Any]) -> dict[str, Any]:
    """Compute expected rows without calling the adapter's helper code."""

    a0 = data["A0"]
    a1 = data["A1"]
    s = data["S"]
    sz = data["S_Z"]
    s_entries = [s[i][i] for i in range(SIZE)]
    t_entries = [
        Fraction(1), data["epsilon"], data["epsilon"],
        data["Pstar_squared"] * data["epsilon"], Fraction(1), data["epsilon"],
    ]
    tz_entries = [Fraction(0) for _ in range(SIZE)]
    si = _inv(s_entries)
    t = _diag(t_entries)
    ti = _inv(t_entries)
    tz = _diag(tz_entries)
    b0 = _add(_matmul(_matmul(si, a0), s), _matmul(_matmul(si, a1), sz))
    b1 = _matmul(_matmul(si, a1), s)
    c0 = _add(_matmul(_matmul(ti, b0), t), _matmul(_matmul(ti, b1), tz))
    c1 = _matmul(_matmul(ti, b1), t)
    qres = _add(_add(b0, data["Akin"], -1), data["Achi"], -1)
    p = [[Fraction(0) for _ in range(SIZE)] for _ in range(SIZE)]
    p[E51[0]][E51[1]] = 2 * data["chi"]
    p_z = [[Fraction(0) for _ in range(SIZE)] for _ in range(SIZE)]
    p_z[E51[0]][E51[1]] = 2 * data["chi_Z"]
    w_matrix = _matmul(s, t)
    wz_matrix = _add(_add(_matmul(sz, t), _matmul(s, tz)), w_matrix)
    y = data["Y"]
    y_z = data["Y_Z"]
    return {
        "B0_star": b0,
        "B1_star": b1,
        "B0_extra": c0,
        "B1_extra": c1,
        "Qres": qres,
        "Qres_with_P": _add(qres, p),
        "P": p,
        "P_Z": p_z,
        "T": t,
        "T_Z": tz,
        "W": _vec(w_matrix, y),
        "W_Z": [a + b for a, b in zip(_vec(wz_matrix, y), _vec(w_matrix, y_z))],
        "P_Y": _vec(p, y),
        "P_Y_Z": [a + b for a, b in zip(_vec(p_z, y), _vec(p, y_z))],
        "reverse_A0": _add(
            _matmul(_matmul(s, b0), si),
            _matmul(_matmul(a1, sz), si),
            -1,
        ),
        "reverse_A1": _matmul(_matmul(s, b1), si),
    }


def run_fixture(build: Callable[..., Mapping[str, Any]]) -> dict[str, Any]:
    """Compare every exported synthetic row against independent Fractions."""

    data = synthetic_inputs()
    conditions = {"synthetic_diagonal_S": True, "source_identity_external": True}
    build_data = {key: value for key, value in data.items()
                  if key not in {"epsilon", "Pstar_squared"}}
    result = build(**build_data, conditions=conditions)
    expected = _reference(data)
    checks: dict[str, bool] = {}
    for key in (
        "B0_star", "B1_star", "B0_extra", "B1_extra", "Qres",
        "Qres_with_P", "P", "P_Z", "T", "T_Z",
    ):
        checks[key] = _same(result.get(key), expected[key])
    state = result.get("state", {})
    checks["state_W"] = _same(state.get("W"), expected["W"])
    checks["state_W_Z"] = _same(state.get("W_Z"), expected["W_Z"])
    checks["state_P_Y"] = _same(state.get("P_Y"), expected["P_Y"])
    checks["state_P_Y_Z"] = _same(state.get("P_Y_Z"), expected["P_Y_Z"])
    checks["reverse_A0"] = _same(result.get("reverse_identity", {}).get("A0_from_S"), expected["reverse_A0"])
    checks["reverse_A1"] = _same(result.get("reverse_identity", {}).get("A1_from_S"), expected["reverse_A1"])
    checks["reverse_passed"] = result.get("reverse_identity", {}).get("passed") is True
    checks["nonconstant_S_retained"] = result.get("gauge_scope", {}).get("S_nonconstant_allowed") is True
    checks["explicit_extra_weights_not_inferred"] = result.get("extra_weights") is None
    checks["real_extra_formula_caller_declared"] = result.get("gauge_scope", {}).get(
        "real_natural_component_formula_caller_declared"
    ) is True
    checks["real_extra_formula_not_proven"] = result.get("gauge_scope", {}).get(
        "real_natural_component_formula_verified"
    ) is False
    checks["independent_Pstar_squared_slot"] = (
        result.get("T", [[None] * SIZE for _ in range(SIZE)])[3][3]
        == data["Pstar_squared"] * data["epsilon"]
    )
    checks["real_extra_T_Z_zero"] = all(
        value == 0
        for row in result.get("T_Z", [])
        for value in row
    )
    checks["extra_principal_commutation"] = result.get("gauge_scope", {}).get(
        "extra_principal_commutes", {}
    ).get("Akin") is True and result.get("gauge_scope", {}).get(
        "extra_principal_commutes", {}
    ).get("Achi") is True
    checks["P_applied_only_after_S1_S5"] = result.get("P_status", {}).get("applied") is True
    checks["P_Z_leibniz"] = result.get("ordinary_Z1_leibniz", {}).get(
        "P_product_rule_available"
    ) is True
    checks["conditions_retained"] = result.get("ownership", {}).get("conditions") is conditions
    checks["source_identity_not_claimed"] = result.get("same_source_identity_proved") is False
    checks["arbitrary_matrix_gauge_rejected_by_scope"] = result.get("gauge_scope", {}).get(
        "arbitrary_matrix_gauge_supported"
    ) is False
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "nonconstant_S": True,
        "state_and_chi_nonzero": True,
        "source_identity_proved": False,
        "method": "independent exact Fraction matrix products and ordinary-Z product rows",
    }


def run(build: Any = None, *_unused: Any) -> dict[str, Any]:
    """Alias used by controlled fixture consumers.

    A module/object exposing ``build`` is accepted in addition to the bare
    callable so a consumer need not manufacture an import-time adapter
    wrapper.  The optional second positional argument is deliberately ignored
    for compatibility with older fixture runners that pass an arithmetic
    context; this fixture itself remains exact ``Fraction`` arithmetic.
    """

    if build is None:
        from . import n1_gauge_conjugation_identity_v1 as adapter
        build = adapter.build
    elif not callable(build) and hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be a conjugation callable or module")
    return run_fixture(build)


__all__ = ["synthetic_inputs", "run_fixture", "run"]
