"""Conditional, non-numerical main/gap first-moment algebra.

This is an unexecuted, source-independent draft.  It deliberately does not
import a frame builder, a gap builder, or a checker, and it does not read
``root``.  The only source inputs in scope are abstract placeholders.

The first part records the finite-``mu`` source-moment implication

    A[C]_k = -mu * (W*M_k + ap_k*W*J/mu)

and checks the gap/main normalisation for ``k = 0, 1, 2``.  The companion
row ``A[C]_k + mu*D_k`` is retained verbatim, but is explicitly outside this
proof scope.  The second part derives ordinary-``Z`` product rules for a
``C2`` product ``m = theta*history`` and records endpoint identities for
``R**2*Q``, ``R**2*Q**2``, ``R**2*U**2`` and ``R*Qy``.

Nothing here binds a real frame seam, an actual ``Z`` transport, a numerical
value, an error budget, or an error reduction.  A caller may later use this
abstract receipt as one algebraic premise only after supplying a separate,
source-bound input package.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


SCHEMA = "frame_main_gap_moment_algebra_v1"
ORDERS = (0, 1, 2)


def _sympy():
    """Import the exact symbolic backend only when a caller enters ``recover``."""
    import sympy as sp

    return sp


def _reduce(sp: Any, expression: Any) -> Any:
    """Apply exact, non-numerical simplification to a symbolic expression."""
    # ``powsimp(..., force=True)`` is needed for the positive ``mu`` scale
    # identity; ``together``/``cancel`` are useful for the quotient endpoint
    # formulas.  No floating-point conversion is used.
    return sp.simplify(sp.cancel(sp.together(sp.powsimp(expression, force=True))))


def _zero(sp: Any, expression: Any) -> bool:
    """Return the exact SymPy zero test used by every proof check."""
    return _reduce(sp, expression) == 0


def _text(expression: Any) -> str:
    """Keep returned evidence JSON-friendly while preserving exact syntax."""
    return str(expression)


def recover(root: Path | str | None) -> dict[str, Any]:
    """Return an abstract exact-algebra receipt without consuming ``root``.

    ``root`` is accepted only to match the repository recovery interface.  It
    is intentionally unused: no archive, source function, saved receipt,
    builder, or checker is consulted by this proof.  SymPy construction is
    lazy so importing this draft itself performs no symbolic evaluation.
    """
    del root
    sp = _sympy()

    # ------------------------------------------------------------------
    # A. Finite-mu source moments and the first-row normalisation.
    # ------------------------------------------------------------------
    mu = sp.Symbol("mu", positive=True)
    lambda_ = sp.Symbol("lambda", real=True)
    logpref = sp.Symbol("logpref", real=True)
    J = sp.Symbol("J", real=True)
    M_0, M_1, M_2 = sp.symbols("M_0 M_1 M_2", real=True)
    ap_0, ap_1, ap_2 = sp.symbols("ap_0 ap_1 ap_2", real=True)
    C_0, C_1, C_2 = sp.symbols("C_0 C_1 C_2", real=True)
    D_0, D_1, D_2 = sp.symbols("D_0 D_1 D_2", real=True)

    M = (M_0, M_1, M_2)
    ap = (ap_0, ap_1, ap_2)
    C = (C_0, C_1, C_2)
    D = (D_0, D_1, D_2)

    W = sp.exp(-13 * lambda_ / mu - logpref)
    H = sp.exp(logpref - sp.log(mu) + 2 * lambda_ / mu)
    E = sp.exp(-11 * lambda_ / mu)

    # A[C]_k is kept as an indexed symbolic source operator.  The C_k labels
    # remain in the returned source schema; no source value is selected.
    A_of_C = sp.IndexedBase("A_C")
    A_C = tuple(A_of_C[k] for k in ORDERS)
    finite_mu_rhs = tuple(
        -mu * (W * M[k] + ap[k] * W * J / mu) for k in ORDERS
    )
    finite_mu_rows = tuple(
        sp.Eq(A_C[k], finite_mu_rhs[k]) for k in ORDERS
    )

    # Preserve the second row instead of silently replacing it by the first.
    # It is deliberately not normalised or checked below.
    second_rows_A_plus_mu_D = tuple(
        A_C[k] + mu * D[k] for k in ORDERS
    )

    gap_normalized_rows = tuple(-H * A_C[k] for k in ORDERS)
    main_normalized_rows = tuple(
        E * (M[k] + ap[k] * J / mu) for k in ORDERS
    )
    gap_from_assumed_source = tuple(-H * finite_mu_rhs[k] for k in ORDERS)
    first_row_residuals = tuple(
        _reduce(sp, gap_from_assumed_source[k] - main_normalized_rows[k])
        for k in ORDERS
    )
    scale_residual = _reduce(sp, H * mu * W - E)

    first_row_checks = {
        "H_mu_W_equals_E": scale_residual == 0,
        **{
            f"first_row_normalization_k{k}_zero": first_row_residuals[k] == 0
            for k in ORDERS
        },
    }

    # ------------------------------------------------------------------
    # B. Complete-Z product jets and the actual Q/Qy quotient algebra.
    # ------------------------------------------------------------------
    Z = sp.Symbol("Z", real=True)
    delta = sp.Symbol("delta", real=True)
    h = 1 - delta
    R = sp.Symbol("R", real=True)

    theta = sp.Function("theta")(Z)
    history = sp.Function("history")(Z)
    m_product = theta * history
    m_product_Z = sp.diff(m_product, Z)
    m_product_ZZ = sp.diff(m_product, Z, 2)
    m_y = -m_product
    U = sp.Integer(0)
    U_y = sp.Integer(0)

    L = 1 - delta * Z**2
    L_Z = sp.diff(L, Z)
    Q_actual = (-h * Z * m_product - (1 - Z**2) * m_product_Z) / L
    # Apply the same linear radial recovery to (Uy, my), where
    # Uy=0 and my=-m. This derives Qy=-Q, rather than assuming it.
    Qy_actual = (-h * Z * m_y - (1 - Z**2) * sp.diff(m_y, Z)) / L

    # A separate generic m(Z) keeps the quotient rule legible and makes the
    # m_ZZ and L_Z terms visible in the explicit Q_Z identity.
    m_function = sp.Function("m")(Z)
    m_Z = sp.diff(m_function, Z)
    m_ZZ = sp.diff(m_function, Z, 2)
    numerator = -h * Z * m_function - (1 - Z**2) * m_Z
    numerator_Z = (
        -h * m_function
        - h * Z * m_Z
        + 2 * Z * m_Z
        - (1 - Z**2) * m_ZZ
    )
    Q = numerator / L
    Qy = -Q
    # Keep a named L_Z atom for the displayed identity, then substitute its
    # actual ordinary derivative for the executable equality check.
    L_Z_atom = sp.Symbol("L_Z")
    Q_Z_explicit_named = (numerator_Z * L - numerator * L_Z_atom) / L**2
    Q_Z_explicit = Q_Z_explicit_named.subs(L_Z_atom, L_Z)

    # Ordinary-Z Leibniz rules are recorded independently for each moment
    # history, since theta is allowed to depend on Z.
    history_by_k = tuple(
        sp.Function(f"history_{k}")(Z) for k in ORDERS
    )
    moment_products = tuple(theta * history_by_k[k] for k in ORDERS)
    moment_leibniz = {}
    for k in ORDERS:
        history_k = history_by_k[k]
        m_k = moment_products[k]
        m_k_Z_expected = sp.diff(theta, Z) * history_k + theta * sp.diff(history_k, Z)
        m_k_ZZ_expected = (
            sp.diff(theta, Z, 2) * history_k
            + 2 * sp.diff(theta, Z) * sp.diff(history_k, Z)
            + theta * sp.diff(history_k, Z, 2)
        )
        moment_leibniz[k] = {
            "m_k": _text(m_k),
            "m_k_Z": _text(sp.diff(m_k, Z)),
            "m_k_Z_expected": _text(m_k_Z_expected),
            "m_k_ZZ": _text(sp.diff(m_k, Z, 2)),
            "m_k_ZZ_expected": _text(m_k_ZZ_expected),
            "m_k_Z_check": _zero(sp, sp.diff(m_k, Z) - m_k_Z_expected),
            "m_k_ZZ_check": _zero(sp, sp.diff(m_k, Z, 2) - m_k_ZZ_expected),
        }

    q_checks = {
        "m_y_equals_minus_m": _zero(sp, m_y + m_product),
        "U_equals_zero": U == 0,
        "U_y_equals_zero": U_y == 0,
        "Qy_equals_minus_Q": _zero(sp, Qy_actual + Q_actual),
        "Q_Z_quotient_identity": _zero(sp, sp.diff(Q, Z) - Q_Z_explicit),
        "Q_Z_explicitly_contains_m_ZZ": Q_Z_explicit_named.has(m_ZZ),
        "Q_Z_explicitly_contains_L_Z": Q_Z_explicit_named.has(L_Z_atom),
    }

    # ------------------------------------------------------------------
    # C. Ordinary Z0/Z1 jets at arbitrary Z, never Z=0 or Z=1 samples.
    # ------------------------------------------------------------------
    m_0, m_1, m_2 = sp.symbols("m_0 m_1 m_2", real=True)
    N0 = -h * Z * m_0 - (1 - Z**2) * m_1
    N1 = -h * (m_0 + Z*m_1) + 2*Z*m_1 - (1 - Z**2)*m_2
    q0 = N0 / L
    q1 = N1 / L - N0 * L_Z / L**2
    R_Z = sp.Integer(0)
    U_Z = sp.Integer(0)
    product_generic = {"R2Q": R**2 * Q, "R2Q2": R**2 * Q**2,
                       "R2U2": R**2 * U**2, "RQy": R * Qy}
    product_Z_generic = {"R2Q": R**2 * Q_Z_explicit,
        "R2Q2": 2*R**2*Q*Q_Z_explicit, "R2U2": sp.Integer(0),
        "RQy": -R*Q_Z_explicit}
    product_leibniz_checks = {name: _zero(sp, sp.diff(product_generic[name], Z) - product_Z_generic[name]) for name in product_generic}
    endpoint_value_expressions = {"R2Q": R**2*q0, "R2Q2": R**2*q0**2,
        "R2U2": sp.Integer(0), "RQy": -R*q0}
    endpoint_Z_expressions = {"R2Q": R**2*q1, "R2Q2": 2*R**2*q0*q1,
        "R2U2": sp.Integer(0), "RQy": -R*q1}
    endpoint_values = {name: _text(value) for name, value in endpoint_value_expressions.items()}
    endpoint_Z_values = {name: _text(value) for name, value in endpoint_Z_expressions.items()}
    endpoint_raw_jets = {"Q": _text(q0), "Q_Z": _text(q1), "Qy": _text(-q0),
        "L": _text(L), "L_Z": _text(L_Z), "coordinate": "arbitrary Z in L!=0 domain"}
    m0_left, m1_left, m2_left = sp.symbols("m_0_left m_1_left m_2_left", real=True)
    m0_right, m1_right, m2_right = sp.symbols("m_0_right m_1_right m_2_right", real=True)
    R_left, R_right = sp.symbols("R_left R_right", real=True)
    left_substitution = {m_0: m0_left, m_1: m1_left, m_2: m2_left, R: R_left}
    right_substitution = {m_0: m0_right, m_1: m1_right, m_2: m2_right, R: R_right}
    same_side_substitution = {m0_right: m0_left, m1_right: m1_left,
        m2_right: m2_left, R_right: R_left}
    side_equality_checks: dict[str, bool] = {}
    for product_name in ("R2Q", "R2Q2", "R2U2", "RQy"):
        for order, expression in (("Z0", endpoint_value_expressions[product_name]),
                                  ("Z1", endpoint_Z_expressions[product_name])):
            left = expression.subs(left_substitution)
            right = expression.subs(right_substitution).subs(same_side_substitution)
            side_equality_checks[f"{product_name}_{order}_same_side"] = _zero(sp, right-left)

    proof_checks = {
        **first_row_checks,
        **{f"moment_product_k{k}_{name}": value
           for k, checks in moment_leibniz.items()
           for name, value in checks.items()
           if name.endswith("_check")},
        **q_checks,
        **{f"product_leibniz_{name}": value
           for name, value in product_leibniz_checks.items()},
        **side_equality_checks,
    }

    return {
        "schema": SCHEMA,
        "status": "conditional_symbolic_abstract_only",
        "passed": all(proof_checks.values()),
        "symbolic_backend": "sympy_exact_forms_constructed_only_when_recover_is_called",
        "proof_checks": proof_checks,
        "assumptions": {
            "mu_positive": True,
            "gp_and_parameters_Z_independent": True,
            "gp_endpoint_gp_11_equals_zero": True,
            "gp_endpoint_gp_prime_11_equals_zero": True,
            "theta_and_history_same_source_C2": True,
            "moment_identity_is_a_complete_Z_function": True,
            "L_nonzero_on_scoped_Z_domain": True,
            "R_Z_equals_zero": True,
            "same_side_m_0_m_1_m_2": True,
            "same_side_R": True,
            "ordinary_Z_rows": "Z0 is the value function and Z1 is its first derivative, both at arbitrary Z",
            "no_axial_sample_values_used": True,
        },
        "finite_mu_algebra": {
            "symbols": {
                "mu": _text(mu),
                "lambda": _text(lambda_),
                "logpref": _text(logpref),
                "M": [_text(value) for value in M],
                "ap": [_text(value) for value in ap],
                "J": _text(J),
                "C": [_text(value) for value in C],
                "D": [_text(value) for value in D],
            },
            "index_meaning": "k=0,1,2 are ordinary Z derivative orders of the first moment row, not distinct radial moment channels",
            "J_interpretation": "formal integral from 0 to 11 of exp(lambda*xi/mu)*gp(xi) dxi; no integral is evaluated or source binding inferred",
            "W": _text(W),
            "H": _text(H),
            "E": _text(E),
            "scale_identity": "H*mu*W = E",
            "scale_identity_residual": _text(scale_residual),
            "source_moment_rows": [_text(row) for row in finite_mu_rows],
            "gap_normalized_rows": [_text(row) for row in gap_normalized_rows],
            "main_normalized_rows": [_text(row) for row in main_normalized_rows],
            "first_row_residuals": [_text(value) for value in first_row_residuals],
            "second_rows_A_plus_mu_D": [_text(row) for row in second_rows_A_plus_mu_D],
            "proof_scope": "first finite-mu row only",
            "second_row_in_scope": False,
            "second_row_preserved": True,
        },
        "ordinary_Z_product_algebra": {
            "m_equals_theta_times_history": _text(m_product),
            "m_Z": _text(m_product_Z),
            "m_ZZ": _text(m_product_ZZ),
            "my": _text(m_y),
            "U": _text(U),
            "Uy": _text(U_y),
            "L": _text(L),
            "L_Z": _text(L_Z),
            "Q": _text(Q_actual),
            "Qy": _text(Qy_actual),
            "Q_generic": _text(Q),
            "Q_Z_explicit_named": _text(Q_Z_explicit_named),
            "Q_Z_explicit": _text(Q_Z_explicit),
            "moment_leibniz": moment_leibniz,
            "product_leibniz_forms": {
                name: _text(expression) for name, expression in product_Z_generic.items()
            },
        },
        "endpoint_forms": {
            "values": endpoint_values,
            "ordinary_Z1": endpoint_Z_values,
            "raw_jets": endpoint_raw_jets,
            "side_equality_scope": {
                "same_m_0": "m(Z) agrees on both sides",
                "same_m_1": "m_Z(Z) agrees on both sides",
                "same_m_2": "m_ZZ(Z) agrees on both sides",
                "R_Z": "0",
                "actual_source_function_identity_claimed": False,
            },
        },
        "source_scope": {
            "root_consumed": False,
            "abstract_input_only": True,
            "old_builder_called": False,
            "new_builder_called": False,
            "old_checker_called": False,
            "new_checker_called": False,
            "real_frame_seam_proved": False,
            "actual_Z_transport_proved": False,
            "numerical_values_selected": False,
            "error_reduction_proved": False,
            "single_run_or_freeze_receipt": False,
            "caller_must_supply_source_bound_package": True,
        },
        "limitations": (
            "Conditional symbolic algebra only. The finite-mu source rows, "
            "same-side jets, C2 product premise, gp endpoint premise, and "
            "L nonvanishing premise are assumptions; no actual frame seam, "
            "physical Z transport, numerical bound, or error reduction is "
            "proved here."
        ),
    }

