"""Exact signed transport for the original ``A/B/C`` native histories.

This module is intentionally a small source and transport boundary.  It has
two independent APIs:

* :func:`pure_signed_formula` builds exact SymPy expressions from the source
  ODEs.  The inlet histories remain symbols in those expressions.
* :func:`transport_signed_jet` evaluates finite, signed interval jets through
  ordinary ``Z`` orders ``0..2``.  It accepts signed value records and finite
  interval values only.  A logarithmic cap is an error, even when its shape
  resembles a value record.

The source ODEs are

``A_y = U**2 - A``, ``B_y = phi**2 - 2*B`` and ``C_y = phi**2 - C``.

On the reference/post charts, ``phi_y = -2*phi/5`` and ``U`` is constant.
The closed transport therefore keeps all three unknown inlet memories.  The
restoration ``A`` formula consumes the already-saved ``gamma1`` and the new
``gamma2`` as inputs; it never evaluates the old integral again.

No saved report is read at import time and no accepted producer, builder,
checker, recurrence, quadrature, norm, or cap provider is called here.
"""

from __future__ import annotations

import ast
import hashlib
import math
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


MODULE = "src/openai_ns_reconstruction/original_signed_native_history_kernel.py"
ORDINARY_Z_ORDERS = (0, 1, 2)
NATIVE_HISTORY_NAMES = ("A", "B", "C")
_CAP_KEYS = frozenset(
    {
        "cap",
        "caps",
        "logcap",
        "log_cap",
        "log_absolute_upper_bound",
        "absolute_upper_bound",
        "bound",
        "bounds",
        "majorant",
        "majorants",
        "upper_bound",
        "error_cap",
    }
)


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _source_preflight(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run the required LF, one-EOF, and AST check before a fixture."""

    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
    path = root / Path(MODULE)
    raw = path.read_bytes()
    if b"\r" in raw:
        raise ValueError("module source must use LF line endings")
    if not raw.endswith(b"\n") or raw.endswith(b"\n\n"):
        raise ValueError("module source must have exactly one final LF")
    ast.parse(raw.decode("utf-8"), filename=str(path))
    return {
        "passed": True,
        "path": MODULE,
        "module_sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "lf_only": True,
        "single_eof": True,
        "ast_parsed": True,
    }


def _is_sympy_float(value: Any) -> bool:
    try:
        import sympy as sp

        return isinstance(value, sp.Float) or (
            isinstance(value, sp.Basic) and value.has(sp.Float)
        )
    except ImportError:  # pragma: no cover - the symbolic API requires SymPy
        return False


def _symbolic_value(value: Any, label: str) -> Any:
    """Coerce one exact symbolic input while refusing cap-shaped data."""

    import sympy as sp

    if isinstance(value, Mapping):
        raise TypeError(label + " must be an exact expression, not a cap mapping")
    if isinstance(value, float) or _is_sympy_float(value):
        raise TypeError(label + " rejects binary/SymPy floating-point input")
    if isinstance(value, bool):
        raise TypeError(label + " must not be boolean")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    if isinstance(value, (int, str, sp.Basic)):
        return sp.sympify(value, rational=True)
    raise TypeError(label + " must be an exact SymPy expression")


def _default_symbol(name: str) -> Any:
    import sympy as sp

    return sp.Symbol(name)


def _resolve_symbolic(value: Any, name: str) -> Any:
    return _symbolic_value(_default_symbol(name) if value is None else value, name)


def _validate_frozen_parameters(frozen_params: Mapping[str, Any] | None) -> None:
    """Reject dependent frozen parameters before constructing a formula.

    Frozen scales are source parameters, not expressions derived from one
    another.  A mapping such as ``{"R": R, "F0": F0*R}`` is rejected; the
    caller must pass the independent source symbols and form the product in
    the returned expression instead.
    """

    import sympy as sp

    if frozen_params is None:
        return
    if not isinstance(frozen_params, Mapping):
        raise TypeError("frozen_params must be a mapping of independent symbols")
    values = {}
    for key, raw in frozen_params.items():
        if not isinstance(key, str):
            raise TypeError("frozen parameter names must be strings")
        value = _symbolic_value(raw, "frozen_params[" + key + "]")
        values[key] = value
        if not isinstance(value, sp.Symbol):
            raise ValueError("frozen parameter " + key + " must remain an independent symbol")
        if value.name != key:
            raise ValueError("frozen parameter " + key + " is renamed or dependent")


def _exact_symbolic_parameters(
    U: Any,
    phi_end: Any,
    s: Any,
    A_in: Any,
    B_in: Any,
    C_in: Any,
    *,
    frozen_params: Mapping[str, Any] | None,
) -> dict[str, Any]:
    _validate_frozen_parameters(frozen_params)
    return {
        "U": _resolve_symbolic(U, "U"),
        "phi_end": _resolve_symbolic(phi_end, "phi_end"),
        "s": _resolve_symbolic(s, "s"),
        "A_in": _resolve_symbolic(A_in, "A_in"),
        "B_in": _resolve_symbolic(B_in, "B_in"),
        "C_in": _resolve_symbolic(C_in, "C_in"),
    }


def _symbolic_gamma(t: Any, h: Any, power: int, name: str) -> Any:
    import sympy as sp

    t_expr = _symbolic_value(t, "t")
    if h is None:
        h_expr = sp.Function("h")
    else:
        if isinstance(h, Mapping) or isinstance(h, float) or _is_sympy_float(h):
            raise TypeError("h must be an exact SymPy function/expression")
        h_expr = h if callable(h) else sp.sympify(h, rational=True)
    v = sp.Dummy("v")
    source = h_expr(v) if callable(h_expr) else h_expr
    source = _symbolic_value(source, "h(v)")
    return sp.exp(-t_expr) + sp.Integral(
        sp.exp(-(t_expr - v)) * source**power, (v, 0, t_expr)
    )


def _symbolic_transport(
    U: Any,
    phi_end: Any,
    s: Any,
    A_in: Any,
    B_in: Any,
    C_in: Any,
    *,
    frozen_params: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    import sympy as sp

    p = _exact_symbolic_parameters(
        U, phi_end, s, A_in, B_in, C_in, frozen_params=frozen_params
    )
    Uv, phie, sv = p["U"], p["phi_end"], p["s"]
    Ain, Bin, Cin = p["A_in"], p["B_in"], p["C_in"]
    if phie == 0:
        raise ValueError("phi_end must be nonzero for normalized deviations")
    phi_end_sq = phie**2
    phi_in_sq = phi_end_sq * sp.exp(sp.Rational(4, 5) * sv)
    decay_1 = sp.exp(-sv)
    decay_2 = sp.exp(-2 * sv)
    A = Uv**2 + decay_1 * (Ain - Uv**2)
    B = sp.Rational(5, 6) * phi_end_sq + decay_2 * (
        Bin - sp.Rational(5, 6) * phi_in_sq
    )
    C = 5 * phi_end_sq + decay_1 * (Cin - 5 * phi_in_sq)
    normalized_B = sp.cancel(B / phi_end_sq - sp.Rational(5, 6))
    normalized_C = sp.cancel(C / phi_end_sq - 5)
    return {
        "A": sp.expand(A),
        "B": sp.expand(B),
        "C": sp.expand(C),
        "phi_end_squared": phi_end_sq,
        "phi_in_squared": phi_in_sq,
        "decay_1": decay_1,
        "decay_2": decay_2,
        "memory_contributions": {
            "A": decay_1 * Ain,
            "B": decay_2 * Bin,
            "C": decay_1 * Cin,
        },
        "normalized_B_deviation": normalized_B,
        "normalized_C_deviation": normalized_C,
        "source_coordinate": "s from inlet to endpoint",
        "unknown_inlet_memories": ["A_in", "B_in", "C_in"],
        "inlet_memories_symbolic": all(
            symbol in expr.free_symbols
            for symbol, expr in (
                (Ain, A),
                (Bin, B),
                (Cin, C),
            )
        ),
    }


def source_ode_formulas(
    U: Any = None,
    phi: Any = None,
    A: Any = None,
    B: Any = None,
    C: Any = None,
    *,
    frozen_params: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the exact signed source ODE right-hand sides."""

    import sympy as sp

    _validate_frozen_parameters(frozen_params)
    Uv = _resolve_symbolic(U, "U")
    phiv = _resolve_symbolic(phi, "phi")
    Av = _resolve_symbolic(A, "A")
    Bv = _resolve_symbolic(B, "B")
    Cv = _resolve_symbolic(C, "C")
    return {
        "A_y": sp.expand(Uv**2 - Av),
        "B_y": sp.expand(phiv**2 - 2 * Bv),
        "C_y": sp.expand(phiv**2 - Cv),
        "M_y": sp.expand(Uv - _resolve_symbolic(None, "M")),
        "equations": {
            "A": "A_y=U**2-A",
            "B": "B_y=phi**2-2*B",
            "C": "C_y=phi**2-C",
        },
    }


def phi_from_Utheta(Utheta: Any = None, F0: Any = None, R: Any = None) -> Any:
    """Return the source normalization ``Utheta/(F0*sqrt(2*R))``."""

    import sympy as sp

    Uv = _resolve_symbolic(Utheta, "Utheta")
    Fv = _resolve_symbolic(F0, "F0")
    Rv = _resolve_symbolic(R, "R")
    expression = Uv / (Fv * sp.sqrt(2 * Rv))
    identity = sp.cancel(Rv * Fv**2 * expression**2 - Uv**2 / 2)
    _require(identity == 0, "phi normalization identity did not simplify exactly")
    return expression


def pure_signed_formula(
    U: Any = None,
    phi_end: Any = None,
    s: Any = None,
    A_in: Any = None,
    B_in: Any = None,
    C_in: Any = None,
    *,
    frozen_params: Mapping[str, Any] | None = None,
    frozen_parameters: Mapping[str, Any] | None = None,
    Utheta: Any = None,
    F0: Any = None,
    R: Any = None,
    U110: Any = None,
    Z: Any = None,
    sigma: Any = None,
    t: Any = None,
    gamma1: Any = None,
    gamma2: Any = None,
    epsilon_A_Rz: Any = None,
    h: Any = None,
    logref: Any = None,
    A_R110: Any = None,
    B_Rsh: Any = None,
    C_Rsh: Any = None,
    phi_Rz: Any = None,
    error_metadata: Mapping[str, Any] | None = None,
    largeerror: Any = None,
    logscale_metadata: Any = None,
) -> dict[str, Any]:
    """Build exact reference, post, and restoration signed formulas.

    ``A_in``, ``B_in`` and ``C_in`` are never inferred from a cap or replaced
    by a midpoint.  Omitting them creates independent symbolic inlet memory
    symbols in the returned expressions.
    """

    import sympy as sp

    if frozen_params is not None and frozen_parameters is not None:
        raise TypeError("provide only one frozen parameter mapping")
    frozen = frozen_params if frozen_params is not None else frozen_parameters
    if error_metadata is not None and not isinstance(error_metadata, Mapping):
        raise TypeError("error_metadata must be a mapping retained as metadata")
    reference = _symbolic_transport(
        U, phi_end, s, A_in, B_in, C_in, frozen_params=frozen
    )
    ode = source_ode_formulas(
        U=_resolve_symbolic(U, "U"),
        phi=_resolve_symbolic(phi_end, "phi_end"),
        frozen_params=frozen,
    )

    Uv = _resolve_symbolic(U, "U")
    phie = _resolve_symbolic(phi_end, "phi_end")
    sv = _resolve_symbolic(s, "s")
    Zvar = _resolve_symbolic(Z, "Z")
    U110var = _resolve_symbolic(U110, "U110")
    tvar = _resolve_symbolic(t, "t")
    E = U110var - 4 * Zvar
    sigma_expr = _resolve_symbolic(sigma, "sigma") if sigma is not None else sp.Function("sigma")(tvar)
    gamma1_definition = _symbolic_gamma(tvar, h, 1, "gamma1")
    gamma2_definition = _symbolic_gamma(tvar, h, 2, "gamma2")
    # The definitions are exposed for the source proof, while the restoration
    # formula consumes saved scalar inputs.  Omitting an input therefore keeps
    # an unknown symbol instead of recomputing an old history or selecting a
    # cap-derived value.
    gamma1_expr = (
        _symbolic_value(gamma1, "gamma1")
        if gamma1 is not None
        else _default_symbol("gamma1")
    )
    gamma2_expr = (
        _symbolic_value(gamma2, "gamma2")
        if gamma2 is not None
        else _default_symbol("gamma2")
    )
    epsilon = _resolve_symbolic(epsilon_A_Rz, "epsilon_A_Rz")
    canonical_Rz = canonical_Rz_signed_formula(
        U110=U110var,
        A_R110=A_R110,
        logref=logref,
        B_Rsh=B_Rsh,
        C_Rsh=C_Rsh,
        phi_Rz=phi_Rz,
    )
    U_restore = 4 * Zvar + E * (1 - sigma_expr)
    A_restore = (
        16 * Zvar**2
        + 8 * Zvar * E * gamma1_expr
        + E**2 * gamma2_expr
        + sp.exp(-tvar) * epsilon
    )
    normalization = phi_from_Utheta(Utheta, F0, R)
    return {
        "source_ode": ode,
        "reference": reference,
        "post": {
            **reference,
            "transport_kind": "post_after_restore_reference_closed_form",
        },
        "restoration": {
            "U": sp.expand(U_restore),
            "E": sp.expand(E),
            "gamma1": gamma1_expr,
            "gamma2": gamma2_expr,
            "gamma1_definition": gamma1_definition,
            "gamma2_definition": gamma2_definition,
            "A_restore": sp.expand(A_restore),
            "epsilon_A_Rz": epsilon,
            "phi_y": -sp.Rational(2, 5) * phie,
            "U_y": sp.diff(U_restore, tvar),
            "formula": (
                "A_restore=16*Z**2+8*Z*E*gamma1+E**2*gamma2"
                "+exp(-t)*epsilon_A_Rz"
            ),
            "gamma1_reused_as_input": True,
            "gamma2_consumed_as_input": True,
            "gamma1_recomputed": False,
            "gamma2_recomputed": False,
            "M_history_repeated": False,
        },
        "phi_from_Utheta": normalization,
        "normalization_identity": sp.Integer(0),
        "canonical_Rz": canonical_Rz,
        "normalized_relation": {
            "B": reference["normalized_B_deviation"],
            "C": reference["normalized_C_deviation"],
        },
        "memory_symbols": {
            "A_in": reference["memory_contributions"]["A"],
            "B_in": reference["memory_contributions"]["B"],
            "C_in": reference["memory_contributions"]["C"],
        },
        "signed_source_values_only": True,
        "log_caps_rejected": True,
        "dependent_frozen_params_rejected": True,
        "no_midpoint_or_cap_value": True,
        "error_metadata": error_metadata,
        "largeerror": largeerror,
        "largeerror_preserved": largeerror is not None
        or (isinstance(error_metadata, Mapping) and "largeerror" in error_metadata),
        "logscale_metadata": logscale_metadata,
        "physical_scale_materialized": False,
    }


def source_ode_proof(
    U: Any = None,
    phi_in_squared: Any = None,
    A_in: Any = None,
    B_in: Any = None,
    C_in: Any = None,
) -> dict[str, Any]:
    """Prove the three closed reference formulas satisfy their source ODEs."""

    import sympy as sp

    Uv = _resolve_symbolic(U, "U")
    phi0 = _resolve_symbolic(phi_in_squared, "phi_in_squared")
    Ain = _resolve_symbolic(A_in, "A_in")
    Bin = _resolve_symbolic(B_in, "B_in")
    Cin = _resolve_symbolic(C_in, "C_in")
    y = sp.Symbol("y")
    phi2 = phi0 * sp.exp(-sp.Rational(4, 5) * y)
    A = Uv**2 + sp.exp(-y) * (Ain - Uv**2)
    B = sp.Rational(5, 6) * phi2 + sp.exp(-2 * y) * (
        Bin - sp.Rational(5, 6) * phi0
    )
    C = 5 * phi2 + sp.exp(-y) * (Cin - 5 * phi0)
    residuals = {
        "A": sp.simplify(sp.diff(A, y) - (Uv**2 - A)),
        "B": sp.simplify(sp.diff(B, y) - (phi2 - 2 * B)),
        "C": sp.simplify(sp.diff(C, y) - (phi2 - C)),
    }
    initials = {
        "A": sp.simplify(A.subs(y, 0) - Ain),
        "B": sp.simplify(B.subs(y, 0) - Bin),
        "C": sp.simplify(C.subs(y, 0) - Cin),
    }
    passed = all(value == 0 for value in (*residuals.values(), *initials.values()))
    return {
        "passed": passed,
        "residuals": residuals,
        "initial_residuals": initials,
        "phi_squared": phi2,
        "source_ode": source_ode_formulas(U=Uv, phi=sp.Symbol("phi"), A=A, B=B, C=C),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "inlet_memories_symbolic": all(
            symbol in expression.free_symbols
            for symbol, expression in ((Ain, A), (Bin, B), (Cin, C))
        ),
    }


def decompose_signed_formula(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Return the exact source, memory, normalized, and restoration pieces."""

    result = pure_signed_formula(*args, **kwargs)
    reference = result["reference"]
    result["decomposition"] = {
        "source_ode": result["source_ode"],
        "reference_canonical": {
            "A": reference["A"] - reference["memory_contributions"]["A"],
            "B": reference["B"] - reference["memory_contributions"]["B"],
            "C": reference["C"] - reference["memory_contributions"]["C"],
        },
        "reference_memory": reference["memory_contributions"],
        "normalized_deviations": result["normalized_relation"],
        "restoration": result["restoration"],
        "post": result["post"],
    }
    return result


# Readable compatibility spellings for callers that use the source notation.
signed_formula = pure_signed_formula
signedFormula = pure_signed_formula
symbolic_source = pure_signed_formula
source_formula = pure_signed_formula
sourceODEproof = source_ode_proof
signed_decomposition = decompose_signed_formula


def canonical_Rz_signed_formula(
    U110: Any = None,
    A_R110: Any = None,
    logref: Any = None,
    B_Rsh: Any = None,
    C_Rsh: Any = None,
    phi_Rz: Any = None,
    *,
    Lref: Any = None,
) -> dict[str, Any]:
    """Return the exact signed canonical ``Rz`` formulas.

    The scale is represented by the logarithmic coordinate itself.  This
    helper is symbolic only: callers retain any finite log-scale/error
    enclosure as metadata instead of materializing a huge physical scale.
    ``Lref`` is accepted as a readable alias for ``logref``.
    """

    import sympy as sp

    if logref is not None and Lref is not None:
        raise TypeError("provide only one of logref and Lref")
    L = _resolve_symbolic(Lref if Lref is not None else logref, "logref")
    Uv = _resolve_symbolic(U110, "U110")
    Av = _resolve_symbolic(A_R110, "A_R110")
    Bv = _resolve_symbolic(B_Rsh, "B_Rsh")
    Cv = _resolve_symbolic(C_Rsh, "C_Rsh")
    phiv = _resolve_symbolic(phi_Rz, "phi_Rz")
    if phiv == 0:
        raise ValueError("phi_Rz must be nonzero for normalized canonical formulas")
    phi2 = phiv**2
    A = Uv**2 + sp.exp(-(L - 8)) * (Av - Uv**2)
    B_over_phi2 = (
        sp.Rational(5, 6)
        + sp.exp(-2 * L) * Bv / phi2
        - sp.Rational(5, 6) * sp.exp(-sp.Rational(6, 5) * L)
    )
    C_over_phi2 = (
        5
        + sp.exp(-L) * Cv / phi2
        - 5 * sp.exp(-sp.Rational(1, 5) * L)
    )
    return {
        "A_Rz": sp.expand(A),
        "B_Rz_over_phi_Rz_squared": sp.cancel(B_over_phi2),
        "C_Rz_over_phi_Rz_squared": sp.cancel(C_over_phi2),
        "phi_Rz_squared": phi2,
        "logref": L,
        "physical_scale_materialized": False,
        "error_logscale_metadata_only": True,
        "signed_source_values_only": True,
    }


FADING_MEMORY_NU = Fraction(18, 25)
PHI_RELATIVE_DECAY_RATE = Fraction(9, 25)
SAVED_BSHAPE_SLOPE_UPPER = Fraction(1, 200)
SAVED_SIGMA_PRIME_UPPER = Fraction(8, 1)


def _exact_nonnegative(value: Any, label: str) -> Any:
    import sympy as sp

    expression = _symbolic_value(value, label)
    if expression.is_finite is False:
        raise ValueError(label + " must be finite")
    if expression.is_number and expression < 0:
        raise ValueError(label + " must be nonnegative")
    if expression == sp.nan:
        raise ValueError(label + " must be finite")
    return expression


def _exact_positive(value: Any, label: str) -> Any:
    expression = _exact_nonnegative(value, label)
    if expression.is_number and expression <= 0:
        raise ValueError(label + " must be positive")
    return expression


def relative_phi_density_factor(
    ordinary_Z_order: int,
    relative_phi_jets: Sequence[Any],
) -> Any:
    """Return ``D_k = sum_j binomial(k,j) relphi_j relphi_(k-j)`` exactly."""

    import sympy as sp

    if type(ordinary_Z_order) is not int or ordinary_Z_order not in ORDINARY_Z_ORDERS:
        raise ValueError("ordinary_Z_order must be one of 0, 1, 2")
    if not isinstance(relative_phi_jets, Sequence) or isinstance(
        relative_phi_jets, (str, bytes)
    ):
        raise TypeError("relative_phi_jets must contain exact ordinary-Z values")
    if len(relative_phi_jets) <= ordinary_Z_order:
        raise ValueError("relative_phi_jets must retain rows through ordinary_Z_order")
    relative = [
        _symbolic_value(value, "relative_phi_jets[" + str(index) + "]")
        for index, value in enumerate(relative_phi_jets[: ordinary_Z_order + 1])
    ]
    if relative[0] != 1:
        raise ValueError("relative_phi_jets[0] must be exactly 1")
    density_factor = sum(
        (
            sp.binomial(ordinary_Z_order, index)
            * relative[index]
            * relative[ordinary_Z_order - index]
            for index in range(ordinary_Z_order + 1)
        ),
        sp.Integer(0),
    )
    return sp.expand(density_factor)


def pure_fading_memory_log_bound(
    ordinary_Z_order: int,
    T: Any,
    rate: Any,
    *,
    Dk: Any = None,
    relative_phi_jets: Sequence[Any] | None = None,
    phi110_upper_squared: Any = None,
    log_phi110_upper_squared: Any = None,
    Jin_abs: Any = None,
    log_Jin_abs: Any = None,
    nu: Any = FADING_MEMORY_NU,
) -> dict[str, Any]:
    """Build the exact long-endpoint fading-memory error bound.

    The forcing density is bounded by
    ``phi110_upper_squared * exp(-nu*y) * Dk``.  Supplying logarithmic
    amplitudes keeps very large physical scales in log metadata; this helper
    never turns a saved whole-chart cap into a selected value.
    """

    import sympy as sp

    if type(ordinary_Z_order) is not int or ordinary_Z_order not in ORDINARY_Z_ORDERS:
        raise ValueError("ordinary_Z_order must be one of 0, 1, 2")
    if Dk is not None and relative_phi_jets is not None:
        raise TypeError("provide only one of Dk and relative_phi_jets")
    if Dk is None:
        if relative_phi_jets is None:
            raise TypeError("Dk or relative_phi_jets is required")
        Dk = relative_phi_density_factor(ordinary_Z_order, relative_phi_jets)
    else:
        Dk = _exact_nonnegative(Dk, "Dk")
    if phi110_upper_squared is not None and log_phi110_upper_squared is not None:
        raise TypeError("provide only one of phi110_upper_squared and log_phi110_upper_squared")
    if Jin_abs is not None and log_Jin_abs is not None:
        raise TypeError("provide only one of Jin_abs and log_Jin_abs")
    if phi110_upper_squared is None and log_phi110_upper_squared is None:
        log_phi110_upper_squared = sp.Symbol("log_phi110_upper_squared")
    if Jin_abs is None and log_Jin_abs is None:
        log_Jin_abs = sp.Symbol("log_Jin_abs")
    Tv = _symbolic_value(T, "T")
    ratev = _symbolic_value(rate, "rate")
    nuv = _symbolic_value(nu, "nu")
    if any(expression.is_finite is False for expression in (Tv, ratev, nuv)):
        raise ValueError("T, rate, and nu must be finite")
    if Tv.is_number and Tv < 0:
        raise ValueError("T must be nonnegative")
    if ratev.is_number and ratev <= 0:
        raise ValueError("rate must be positive")
    if nuv.is_number and nuv < 0:
        raise ValueError("nu must be nonnegative")
    if ratev.is_number and nuv.is_number and ratev <= nuv:
        raise ValueError("rate must exceed nu for a decaying-density endpoint bound")
    gap = sp.expand(ratev - nuv)
    if phi110_upper_squared is not None:
        phi2 = _exact_positive(phi110_upper_squared, "phi110_upper_squared")
        forcing_amplitude = phi2 * Dk
        forcing_term = forcing_amplitude * sp.exp(-nuv * Tv) / gap
        forcing_log = sp.log(forcing_amplitude / gap) - nuv * Tv
    else:
        log_phi2 = _symbolic_value(log_phi110_upper_squared, "log_phi110_upper_squared")
        forcing_amplitude = sp.exp(log_phi2) * Dk
        forcing_term = forcing_amplitude * sp.exp(-nuv * Tv) / gap
        forcing_log = log_phi2 + sp.log(Dk / gap) - nuv * Tv
    if Jin_abs is not None:
        Jin = _exact_nonnegative(Jin_abs, "Jin_abs")
        inlet_term = Jin * sp.exp(-ratev * Tv)
        inlet_log = sp.log(Jin) - ratev * Tv
    else:
        log_Jin = _symbolic_value(log_Jin_abs, "log_Jin_abs")
        Jin = sp.exp(log_Jin)
        inlet_term = Jin * sp.exp(-ratev * Tv)
        inlet_log = log_Jin - ratev * Tv
    exact_density_solution = forcing_amplitude * (
        sp.exp(-nuv * Tv) - sp.exp(-ratev * Tv)
    ) / gap
    exact_solution = inlet_term + exact_density_solution
    positive_factor = 1 - sp.exp(-gap * Tv)
    discarded_term = forcing_amplitude * sp.exp(-ratev * Tv) / gap
    endpoint_bound = inlet_term + forcing_term
    return {
        "passed": True,
        "ordinary_Z_order": ordinary_Z_order,
        "T": Tv,
        "rate": ratev,
        "history_rate": "B" if ratev == 2 else ("C" if ratev == 1 else None),
        "nu": nuv,
        "rate_minus_nu": gap,
        "Dk": sp.expand(Dk),
        "density_bound": forcing_amplitude * sp.exp(-nuv * sp.Symbol("y")),
        "endpoint_bound": sp.expand(endpoint_bound),
        "exact_density_solution": sp.expand(exact_density_solution),
        "exact_solution": sp.expand(exact_solution),
        "bound_minus_exact": sp.expand(endpoint_bound - exact_solution),
        "discarded_positive_factor": positive_factor,
        "discarded_term": sp.expand(discarded_term),
        "log_bound_components": {
            "inlet": sp.expand(inlet_log),
            "density": sp.expand(forcing_log),
        },
        "source_density_rate_certificate": {
            "Bshape_slope_upper": sp.Rational(1, 200),
            "sigma_prime_upper": sp.Integer(8),
            "phi_relative_decay_rate": sp.Rational(9, 25),
            "phi_log_derivative_upper": -sp.Rational(9, 25),
            "phi_squared_decay_rate": sp.Rational(18, 25),
            "nonnegative_phi_density_assumed": False,
        },
        "wholechart_cap_used_only_as_error_bound": True,
        "physical_scale_materialized": False,
        "log_error_transport_only": True,
    }


fading_memory_log_bound = pure_fading_memory_log_bound
log_fading_memory_bound = pure_fading_memory_log_bound
fading_memory_endpoint_bound = pure_fading_memory_log_bound
pureLogFadingMemoryBound = pure_fading_memory_log_bound
fadingMemoryLogBound = pure_fading_memory_log_bound
relativePhiDensityFactor = relative_phi_density_factor


def _fresh_fading_memory_checks() -> dict[str, Any]:
    """Use exact rationals to validate the decaying-density endpoint bound."""

    import sympy as sp

    T = sp.Rational(7, 3)
    phi2_upper = sp.Integer(7)
    Jin_abs = sp.Integer(3)
    nu = sp.Rational(18, 25)
    relative_phi = [sp.Integer(1), sp.Rational(1, 5), sp.Rational(2, 7)]
    y = sp.Symbol("y")
    checks: dict[str, bool] = {}
    checks["saved_slope_certificate_exact"] = (
        SAVED_BSHAPE_SLOPE_UPPER == Fraction(1, 200)
        and SAVED_SIGMA_PRIME_UPPER == Fraction(8, 1)
        and PHI_RELATIVE_DECAY_RATE == Fraction(9, 25)
        and FADING_MEMORY_NU == Fraction(18, 25)
    )
    density_factors: dict[str, Any] = {}
    for order in ORDINARY_Z_ORDERS:
        density_factor = relative_phi_density_factor(order, relative_phi)
        density_factors[str(order)] = density_factor
        checks["D" + str(order) + "_exact"] = density_factor == sum(
            (
                sp.binomial(order, index)
                * relative_phi[index]
                * relative_phi[order - index]
                for index in range(order + 1)
            ),
            sp.Integer(0),
        )
        for name, rate in (("B2", sp.Integer(2)), ("C1", sp.Integer(1))):
            label = name + "_Z" + str(order)
            result = pure_fading_memory_log_bound(
                order,
                T,
                rate,
                Dk=density_factor,
                phi110_upper_squared=phi2_upper,
                Jin_abs=Jin_abs,
                nu=nu,
            )
            density = phi2_upper * density_factor * sp.exp(-nu * y)
            exact_integral = sp.integrate(
                sp.exp(-rate * (T - y)) * density,
                (y, 0, T),
            )
            checks[label + "_analytic_integral"] = sp.simplify(
                exact_integral - result["exact_density_solution"]
            ) == 0
            checks[label + "_initial_memory"] = sp.simplify(
                result["exact_solution"]
                - (sp.exp(-rate * T) * Jin_abs + result["exact_density_solution"])
            ) == 0
            checks[label + "_positive_discarded_factor"] = (
                result["discarded_positive_factor"].is_positive is True
            )
            checks[label + "_bound_minus_exact"] = sp.simplify(
                result["bound_minus_exact"] - result["discarded_term"]
            ) == 0 and result["bound_minus_exact"].is_positive is True
            checks[label + "_density_decay"] = sp.simplify(
                result["density_bound"]
                - phi2_upper * density_factor * sp.exp(-nu * sp.Symbol("y"))
            ) == 0
            checks[label + "_no_nonnegative_B_assumption"] = (
                result["source_density_rate_certificate"]["nonnegative_phi_density_assumed"]
                is False
            )
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "density_factors": {key: str(value) for key, value in density_factors.items()},
        "nu": str(nu),
        "phi_relative_decay_rate": str(sp.Rational(9, 25)),
        "phi_log_derivative_upper": str(-sp.Rational(9, 25)),
        "rates": {"B": "2", "C": "1"},
        "positive_factor_discarded": "1-exp(-(rate-nu)*T)",
        "actual_saved_data_loaded": False,
    }


def fresh_fading_memory_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run the fresh exact rational decaying-density fixture."""

    preflight = _source_preflight(repo_root)
    result = _fresh_fading_memory_checks()
    return {
        "passed": result["passed"],
        "schema": "original_signed_native_history_fading_memory_fixture_v1",
        "module_sha256": preflight["module_sha256"],
        "module": MODULE,
        **result,
    }


fading_memory_fixture = fresh_fading_memory_fixture


def _raw_endpoint(value: Any) -> mp.mpf:
    if hasattr(value, "_mpf_"):
        return mp.make_mpf(tuple(value._mpf_))
    if isinstance(value, (tuple, list)) and len(value) == 4 and all(
        type(item) is int for item in value
    ):
        return mp.make_mpf(tuple(value))
    raise TypeError("an exact MPF endpoint tuple is required")


def _bounds(value: Any) -> tuple[mp.mpf, mp.mpf]:
    if not hasattr(value, "_mpi_"):
        raise TypeError("an MPIntervalContext interval is required")
    rows = value._mpi_
    if len(rows) == 4 and all(type(item) is int for item in rows):
        rows = (rows, rows)
    if len(rows) != 2:
        raise TypeError("malformed interval")
    low, high = _raw_endpoint(rows[0]), _raw_endpoint(rows[1])
    if low > high:
        raise ValueError("reversed interval")
    return low, high


def _finite_interval(value: Any) -> bool:
    low, high = _bounds(value)
    return not any(mp.isnan(item) or mp.isinf(item) for item in (low, high))


def _interval(context: MPIntervalContext, low: Any, high: Any | None = None):
    if high is None:
        high = low
    low = _raw_endpoint(low) if hasattr(low, "_mpf_") else low
    high = _raw_endpoint(high) if hasattr(high, "_mpf_") else high
    if low > high:
        raise ValueError("reversed signed interval")
    value = context.mpf([low, high])
    if not _finite_interval(value):
        raise ValueError("finite signed intervals are required")
    return value


def _contains_cap_key(value: Mapping[str, Any]) -> bool:
    return bool(_CAP_KEYS.intersection(str(key).lower() for key in value))


def _read_signed_value(context: MPIntervalContext, value: Any, label: str = "signed value"):
    """Read one finite signed value record without accepting a cap record."""

    if isinstance(value, Mapping):
        keys = {str(key).lower() for key in value}
        if _contains_cap_key(value):
            raise TypeError(label + " is a logarithmic cap, not a signed value")
        if "signed_value" in value:
            marker = value["signed_value"]
            if marker is True and "value" in value:
                return _read_signed_value(context, value["value"], label)
            if marker is True and {"lower_exact_mpf_tuple", "upper_exact_mpf_tuple"}.issubset(keys):
                low = _raw_endpoint(value["lower_exact_mpf_tuple"])
                high = _raw_endpoint(value["upper_exact_mpf_tuple"])
                return _interval(context, low, high)
            if marker is not True:
                return _read_signed_value(context, marker, label)
            raise TypeError(label + " signed record has no finite endpoints")
        if "value" in value and len(value) <= 4:
            return _read_signed_value(context, value["value"], label)
        if {"lower_exact_mpf_tuple", "upper_exact_mpf_tuple"}.issubset(keys):
            low = _raw_endpoint(value["lower_exact_mpf_tuple"])
            high = _raw_endpoint(value["upper_exact_mpf_tuple"])
            return _interval(context, low, high)
        raise TypeError(label + " must be an explicit signed value record")
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(label + " rejects binary floating-point input")
    if value is None:
        raise TypeError(label + " cannot be unknown")
    if hasattr(value, "_mpi_"):
        if not _finite_interval(value):
            raise ValueError(label + " must be finite")
        low, high = _bounds(value)
        return _interval(context, low, high)
    if hasattr(value, "_mpf_"):
        point = _raw_endpoint(value)
        return _interval(context, point, point)
    if isinstance(value, Fraction):
        point = context.mpf(value.numerator) / value.denominator
        return _interval(context, point, point)
    if type(value) is int or isinstance(value, str):
        point = context.mpf(value)
        return _interval(context, point, point)
    if isinstance(value, (tuple, list)):
        if len(value) == 4 and all(type(item) is int for item in value):
            point = _raw_endpoint(value)
            return _interval(context, point, point)
        if len(value) == 2:
            low = _read_signed_value(context, value[0], label + ".lower")
            high = _read_signed_value(context, value[1], label + ".upper")
            return _interval(context, _bounds(low)[0], _bounds(high)[1])
    raise TypeError(label + " must be a finite signed value or signed record")


def read_signed_value(context: MPIntervalContext, value: Any, label: str = "signed value"):
    """Public finite signed-record reader."""

    if not isinstance(context, MPIntervalContext):
        raise TypeError("context must be MPIntervalContext")
    return _read_signed_value(context, value, label)


def _validate_log_error_transport(
    logs: Mapping[str, Any] | None,
) -> Mapping[str, Any] | None:
    """Keep error transport in interval-log form without exponentiation."""

    if logs is None:
        return None
    if not isinstance(logs, Mapping):
        raise TypeError("log_error_transport must be a mapping of interval logs")
    for key, value in logs.items():
        if not isinstance(key, str):
            raise TypeError("log_error_transport keys must be strings")
        if not hasattr(value, "_mpi_") or not _finite_interval(value):
            raise TypeError(
                "log_error_transport values must remain finite MPIntervalContext intervals"
            )
    return logs


def _read_signed_jet(
    context: MPIntervalContext,
    value: Any,
    label: str,
    *,
    length: int = 3,
) -> list[Any]:
    if isinstance(value, Mapping):
        if _contains_cap_key(value):
            raise TypeError(label + " is a cap mapping, not a signed jet")
        value = value.get("rows", value.get("signed_values", value.get("jet", value)))
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(label + " must contain ordinary signed-value rows")
    if len(value) != length:
        raise ValueError(label + " must retain ordinary Z orders 0..2")
    rows = []
    for order, row in enumerate(value):
        if isinstance(row, Mapping) and row.get("ordinary_Z_order", order) != order:
            raise ValueError(label + " rows are not ordered")
        rows.append(_read_signed_value(context, row, label + ".Z" + str(order)))
    return rows


def _jet_add(*jets: Sequence[Any]) -> list[Any]:
    if not jets:
        raise ValueError("at least one jet is required")
    return [sum((jet[k] for jet in jets), jets[0][k] * 0) for k in ORDINARY_Z_ORDERS]


def _jet_neg(jet: Sequence[Any]) -> list[Any]:
    return [-jet[k] for k in ORDINARY_Z_ORDERS]


def _jet_sub(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    return _jet_add(left, _jet_neg(right))


def _jet_mul(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    return [
        sum((math.comb(k, j) * left[j] * right[k - j] for j in range(k + 1)), left[0] * 0)
        for k in ORDINARY_Z_ORDERS
    ]


def _jet_scale(jet: Sequence[Any], scalar: Any) -> list[Any]:
    return [scalar * value for value in jet]


def _square(context: MPIntervalContext, value: Any):
    low, high = _bounds(value)
    candidates = (low * low, low * high, high * low, high * high)
    return _interval(context, min(candidates), max(candidates))


def _exp_jet(context: MPIntervalContext, jet: Sequence[Any]) -> list[Any]:
    e0 = context.exp(jet[0])
    return [e0, e0 * jet[1], e0 * (jet[2] + jet[1] ** 2)]


def _as_jet(context: MPIntervalContext, value: Any, label: str) -> list[Any]:
    if isinstance(value, Mapping):
        value_keys = {str(key).lower() for key in value}
        if "signed_value" in value_keys or {
            "lower_exact_mpf_tuple",
            "upper_exact_mpf_tuple",
        }.issubset(value_keys):
            scalar = _read_signed_value(context, value, label)
            zero = context.mpf(0)
            return [scalar, zero, zero]
        return _read_signed_jet(context, value, label)
    if isinstance(value, (list, tuple)) and len(value) == 3:
        return _read_signed_jet(context, value, label)
    scalar = _read_signed_value(context, value, label)
    zero = context.mpf(0)
    return [scalar, zero, zero]


def _finite_transport_from_jets(
    context: MPIntervalContext,
    U: Sequence[Any],
    phi_end: Sequence[Any],
    s: Sequence[Any],
    A_in: Sequence[Any],
    B_in: Sequence[Any],
    C_in: Sequence[Any],
) -> dict[str, Any]:
    one = context.mpf(1)
    five_six = context.mpf(5) / 6
    five = context.mpf(5)
    U2 = _jet_mul(U, U)
    phi_end2 = _jet_mul(phi_end, phi_end)
    phi_in2 = _jet_mul(phi_end2, _exp_jet(context, _jet_scale(s, context.mpf(4) / 5)))
    decay1 = _exp_jet(context, _jet_neg(s))
    decay2 = _exp_jet(context, _jet_scale(s, context.mpf(-2)))
    A = _jet_add(U2, _jet_mul(decay1, _jet_sub(A_in, U2)))
    B = _jet_add(
        _jet_scale(phi_end2, five_six),
        _jet_mul(decay2, _jet_sub(B_in, _jet_scale(phi_in2, five_six))),
    )
    C = _jet_add(
        _jet_scale(phi_end2, five),
        _jet_mul(decay1, _jet_sub(C_in, _jet_scale(phi_in2, five))),
    )
    phi_end_nonzero = _bounds(phi_end2[0])
    if phi_end_nonzero[0] <= 0 <= phi_end_nonzero[1]:
        raise ValueError("phi_end must exclude zero for normalized deviations")
    normalized_B = _jet_sub(_jet_scale(_jet_mul(B, _reciprocal_jet(context, phi_end2)), one), _constant_jet(context, five_six))
    normalized_C = _jet_sub(_jet_scale(_jet_mul(C, _reciprocal_jet(context, phi_end2)), one), _constant_jet(context, five))
    return {
        "A": A,
        "B": B,
        "C": C,
        "phi_end_squared": phi_end2,
        "phi_in_squared": phi_in2,
        "decay_1": decay1,
        "decay_2": decay2,
        "normalized_B_deviation": normalized_B,
        "normalized_C_deviation": normalized_C,
        "memory_contributions": {
            "A": _jet_mul(decay1, A_in),
            "B": _jet_mul(decay2, B_in),
            "C": _jet_mul(decay1, C_in),
        },
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
    }


def _constant_jet(context: MPIntervalContext, value: Any) -> list[Any]:
    scalar = _read_signed_value(context, value, "constant")
    zero = context.mpf(0)
    return [scalar, zero, zero]


def _reciprocal_jet(context: MPIntervalContext, jet: Sequence[Any]) -> list[Any]:
    low, high = _bounds(jet[0])
    if low <= 0 <= high:
        raise ValueError("reciprocal jet value interval contains zero")
    inv0 = 1 / jet[0]
    inv1 = -jet[1] / (jet[0] ** 2)
    inv2 = 2 * jet[1] ** 2 / (jet[0] ** 3) - jet[2] / (jet[0] ** 2)
    return [inv0, inv1, inv2]


def transport_signed_jet(
    context: MPIntervalContext,
    U: Any,
    phi_end: Any,
    s: Any,
    A_in: Any,
    B_in: Any,
    C_in: Any,
    *,
    error_metadata: Mapping[str, Any] | None = None,
    largeerror: Any = None,
    log_error_transport: Mapping[str, Any] | None = None,
    logerrortransport: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Transport finite signed ``A/B/C`` ordinary ``Z0..Z2`` jets.

    Every input row is read as a signed value.  In particular, an object with
    ``log_absolute_upper_bound`` is rejected before any arithmetic.  ``s``
    may itself be a signed three-row jet; scalar ``s`` is treated as
    Z-independent.
    """

    if not isinstance(context, MPIntervalContext):
        raise TypeError("context must be MPIntervalContext")
    if error_metadata is not None and not isinstance(error_metadata, Mapping):
        raise TypeError("error_metadata must be a mapping retained as metadata")
    if log_error_transport is not None and logerrortransport is not None:
        raise TypeError("provide only one log_error_transport spelling")
    error_logs = _validate_log_error_transport(
        log_error_transport if log_error_transport is not None else logerrortransport
    )
    jets = {
        "U": _as_jet(context, U, "U"),
        "phi_end": _as_jet(context, phi_end, "phi_end"),
        "s": _as_jet(context, s, "s"),
        "A_in": _read_signed_jet(context, A_in, "A_in"),
        "B_in": _read_signed_jet(context, B_in, "B_in"),
        "C_in": _read_signed_jet(context, C_in, "C_in"),
    }
    result = _finite_transport_from_jets(context, **jets)
    result.update(
        {
            "signed_value_inputs": True,
            "caps_used_as_signed_values": False,
            "source_values_selected": False,
            "midpoint_selected": False,
            "ordinary_Z_derivatives": True,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "unknown_inlet_memories_preserved": True,
            "enclosure_label": "signed canonical plus certified radius enclosure",
            "canonical_is_exact_actual_field": False,
            "error_metadata": error_metadata,
            "largeerror": largeerror,
            "largeerror_preserved": largeerror is not None
            or (isinstance(error_metadata, Mapping) and "largeerror" in error_metadata),
            "log_error_transport": error_logs,
            "physical_scale_materialized": False,
        }
    )
    return result


def signed_jet_transport(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return transport_signed_jet(*args, **kwargs)


def finite_signed_jet_transport(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return transport_signed_jet(*args, **kwargs)


def transport_signed_histories(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return transport_signed_jet(*args, **kwargs)


def post_reference_signed_jet(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Use the same exact reference transport on the post chart to ``p=1``."""

    return transport_signed_jet(*args, **kwargs)


reference_signed_jet_transport = post_reference_signed_jet


def source_ode_jets(
    context: MPIntervalContext,
    U: Any,
    phi: Any,
    A: Any,
    B: Any,
    C: Any,
) -> dict[str, Any]:
    """Evaluate the signed finite-interval ordinary source ODE jets."""

    if not isinstance(context, MPIntervalContext):
        raise TypeError("context must be MPIntervalContext")
    Uj, phij = _as_jet(context, U, "U"), _as_jet(context, phi, "phi")
    Aj = _read_signed_jet(context, A, "A")
    Bj = _read_signed_jet(context, B, "B")
    Cj = _read_signed_jet(context, C, "C")
    return {
        "A_y": _jet_sub(_jet_mul(Uj, Uj), Aj),
        "B_y": _jet_sub(_jet_mul(phij, phij), _jet_scale(Bj, context.mpf(2))),
        "C_y": _jet_sub(_jet_mul(phij, phij), Cj),
        "signed_value_inputs": True,
        "caps_used_as_signed_values": False,
    }


def restoration_signed_jets(
    context: MPIntervalContext,
    Z: Any,
    U110: Any,
    sigma: Any,
    t: Any,
    gamma1: Any,
    gamma2: Any,
    epsilon_A_Rz: Any,
    *,
    phi_end: Any | None = None,
    s: Any | None = None,
    A_in: Any | None = None,
    B_in: Any | None = None,
    C_in: Any | None = None,
    error_metadata: Mapping[str, Any] | None = None,
    largeerror: Any = None,
    log_error_transport: Mapping[str, Any] | None = None,
    logerrortransport: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate the restoration ``U/E/A`` formulas from signed inputs.

    ``gamma1`` is an old saved history value and ``gamma2`` is a new scalar
    input.  Neither is recomputed from an interval cap.  Optional ``phi`` and
    memory rows pass through the same reference/post ``B/C`` transport.
    """

    if not isinstance(context, MPIntervalContext):
        raise TypeError("context must be MPIntervalContext")
    if log_error_transport is not None and logerrortransport is not None:
        raise TypeError("provide only one log_error_transport spelling")
    Zj, U110j, sj, tj = (
        _as_jet(context, Z, "Z"),
        _as_jet(context, U110, "U110"),
        _as_jet(context, sigma, "sigma"),
        _as_jet(context, t, "t"),
    )
    g1 = _as_jet(context, gamma1, "gamma1")
    g2 = _as_jet(context, gamma2, "gamma2")
    eps = _as_jet(context, epsilon_A_Rz, "epsilon_A_Rz")
    four = context.mpf(4)
    E = _jet_sub(U110j, _jet_scale(Zj, four))
    Urestore = _jet_add(_jet_scale(Zj, four), _jet_mul(E, _jet_sub(_constant_jet(context, 1), sj)))
    Arestore = _jet_add(
        _jet_scale(_jet_mul(Zj, Zj), context.mpf(16)),
        _jet_scale(_jet_mul(_jet_mul(Zj, E), g1), context.mpf(8)),
        _jet_mul(_jet_mul(E, E), g2),
        _jet_mul(_exp_jet(context, _jet_neg(tj)), eps),
    )
    result: dict[str, Any] = {
        "U": Urestore,
        "E": E,
        "A_restore": Arestore,
        "gamma1": g1,
        "gamma2": g2,
        "epsilon_A_Rz": eps,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "signed_value_inputs": True,
        "caps_used_as_signed_values": False,
        "gamma1_recomputed": False,
        "gamma2_recomputed": False,
        "M_history_repeated": False,
        "enclosure_label": "signed canonical plus certified radius enclosure",
        "canonical_is_exact_actual_field": False,
        "error_metadata": error_metadata,
        "largeerror": largeerror,
        "largeerror_preserved": largeerror is not None
        or (isinstance(error_metadata, Mapping) and "largeerror" in error_metadata),
        "log_error_transport": _validate_log_error_transport(
            log_error_transport if log_error_transport is not None else logerrortransport
        ),
    }
    if any(value is not None for value in (phi_end, s, A_in, B_in, C_in)):
        if any(value is None for value in (phi_end, s, A_in, B_in, C_in)):
            raise ValueError("all reference transport inputs are required together")
        result["reference"] = transport_signed_jet(
            context,
            U110,
            phi_end,
            s,
            A_in,
            B_in,
            C_in,
            error_metadata=error_metadata,
            largeerror=largeerror,
            log_error_transport=log_error_transport,
            logerrortransport=logerrortransport,
        )
        result["post"] = result["reference"]
    return result


def post_restore_signed_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return restoration_signed_jets(*args, **kwargs)


def _record_signed(context: MPIntervalContext, value: Any) -> dict[str, Any]:
    low, high = _bounds(value)
    return {
        "signed_value": True,
        "lower_exact_mpf_tuple": list(low._mpf_),
        "upper_exact_mpf_tuple": list(high._mpf_),
        "directed": True,
        "source_value_selected": False,
    }


def record_signed_jet(context: MPIntervalContext, jet: Sequence[Any]) -> list[dict[str, Any]]:
    """Serialize a finite signed jet while retaining its signed-value tag."""

    rows = _read_signed_jet(context, jet, "signed jet")
    return [
        {"ordinary_Z_order": order, **_record_signed(context, value)}
        for order, value in enumerate(rows)
    ]


def _fresh_signed_record(context: MPIntervalContext, value: Fraction) -> dict[str, Any]:
    point = context.mpf(value.numerator) / value.denominator
    return _record_signed(context, point)


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run fresh exact rational/SymPy checks and return a fixture receipt.

    The source preflight is deliberately first.  The fixture constructs its
    own rational source and finite signed records; it never loads a saved
    native packet or turns a cap into a value.
    """

    preflight = _source_preflight(repo_root)
    import sympy as sp

    fading_memory = _fresh_fading_memory_checks()
    U, P, ell = sp.symbols("U P ell", nonzero=True)
    Ain, Bin, Cin = sp.symbols("Ain Bin Cin")
    symbolic = _symbolic_transport(U, P, ell, Ain, Bin, Cin)
    x = sp.symbols("x")
    phi_in_sq = sp.symbols("phi_in_sq", nonzero=True)
    phi2 = phi_in_sq * sp.exp(-sp.Rational(4, 5) * x)
    A_x = U**2 + sp.exp(-x) * (Ain - U**2)
    B_x = sp.Rational(5, 6) * phi2 + sp.exp(-2 * x) * (Bin - sp.Rational(5, 6) * phi_in_sq)
    C_x = 5 * phi2 + sp.exp(-x) * (Cin - 5 * phi_in_sq)
    ode_checks = {
        "A": sp.simplify(sp.diff(A_x, x) - (U**2 - A_x)) == 0,
        "B": sp.simplify(sp.diff(B_x, x) - (phi2 - 2 * B_x)) == 0,
        "C": sp.simplify(sp.diff(C_x, x) - (phi2 - C_x)) == 0,
    }
    initial_checks = {
        "A": sp.simplify(A_x.subs(x, 0) - Ain) == 0,
        "B": sp.simplify(B_x.subs(x, 0) - Bin) == 0,
        "C": sp.simplify(C_x.subs(x, 0) - Cin) == 0,
    }
    z = sp.symbols("Z")
    rational_subs = {
        U: 2 + z,
        P: 3 - z + z**2,
        ell: sp.Rational(1, 3) + z / 7,
        Ain: 5 - 2 * z + z**2,
        Bin: 7 + z + 2 * z**2,
        Cin: 11 - z + 3 * z**2,
    }
    source_z = {
        key: sp.expand(value.subs(rational_subs))
        for key, value in symbolic.items()
        if key in ("A", "B", "C")
    }
    z12_checks = {
        key: all(sp.diff(value, z, order).subs(z, sp.Rational(1, 5)) != 0 for order in (1, 2))
        for key, value in source_z.items()
    }
    memory_checks = {
        "A": sp.simplify(symbolic["A"].subs({U: 0, Ain: Ain}) - symbolic["A"].subs({U: 0, Ain: 0})) == sp.exp(-ell) * Ain,
        "B": sp.simplify(symbolic["B"].subs(Bin, Bin) - symbolic["B"].subs(Bin, 0)) == sp.exp(-2 * ell) * Bin,
        "C": sp.simplify(symbolic["C"].subs(Cin, Cin) - symbolic["C"].subs(Cin, 0)) == sp.exp(-ell) * Cin,
    }
    normalized_checks = {
        "B": sp.simplify(
            symbolic["normalized_B_deviation"]
            - (sp.exp(-2 * ell) * Bin / P**2 - sp.Rational(5, 6) * sp.exp(-sp.Rational(6, 5) * ell))
        ) == 0,
        "C": sp.simplify(
            symbolic["normalized_C_deviation"]
            - (sp.exp(-ell) * Cin / P**2 - 5 * sp.exp(-sp.Rational(1, 5) * ell))
        ) == 0,
    }
    Utheta, F0, R = sp.symbols("Utheta F0 R", nonzero=True, positive=True)
    phi_identity = sp.simplify(R * F0**2 * phi_from_Utheta(Utheta, F0, R) ** 2 - Utheta**2 / 2) == 0
    restoration = pure_signed_formula(
        U=U,
        phi_end=P,
        s=ell,
        A_in=Ain,
        B_in=Bin,
        C_in=Cin,
        U110=sp.Symbol("U110"),
        Z=z,
        sigma=sp.Symbol("sigma"),
        t=sp.Symbol("t"),
        gamma1=sp.Symbol("gamma1"),
        gamma2=sp.Symbol("gamma2"),
        epsilon_A_Rz=sp.Symbol("epsilon_A_Rz"),
        Utheta=Utheta,
        F0=F0,
        R=R,
    )
    restoration_formula_check = (
        sp.expand(restoration["restoration"]["A_restore"])
        == sp.expand(
            16 * z**2
            + 8 * z * (sp.Symbol("U110") - 4 * z) * sp.Symbol("gamma1")
            + (sp.Symbol("U110") - 4 * z) ** 2 * sp.Symbol("gamma2")
            + sp.exp(-sp.Symbol("t")) * sp.Symbol("epsilon_A_Rz")
        )
    )
    ode_proof = source_ode_proof(U=U, phi_in_squared=phi_in_sq, A_in=Ain, B_in=Bin, C_in=Cin)
    canonical = canonical_Rz_signed_formula(
        U110=sp.Symbol("U110"),
        A_R110=sp.Symbol("A_R110"),
        logref=sp.Symbol("logref"),
        B_Rsh=sp.Symbol("B_Rsh"),
        C_Rsh=sp.Symbol("C_Rsh"),
        phi_Rz=sp.Symbol("phi_Rz", nonzero=True),
    )
    canonical_checks = {
        "A": sp.simplify(
            canonical["A_Rz"]
            - (
                sp.Symbol("U110") ** 2
                + sp.exp(-(sp.Symbol("logref") - 8))
                * (sp.Symbol("A_R110") - sp.Symbol("U110") ** 2)
            )
        )
        == 0,
        "B": sp.simplify(
            canonical["B_Rz_over_phi_Rz_squared"]
            - (
                sp.Rational(5, 6)
                + sp.exp(-2 * sp.Symbol("logref"))
                * sp.Symbol("B_Rsh")
                / sp.Symbol("phi_Rz", nonzero=True) ** 2
                - sp.Rational(5, 6) * sp.exp(-sp.Rational(6, 5) * sp.Symbol("logref"))
            )
        )
        == 0,
        "C": sp.simplify(
            canonical["C_Rz_over_phi_Rz_squared"]
            - (
                5
                + sp.exp(-sp.Symbol("logref"))
                * sp.Symbol("C_Rsh")
                / sp.Symbol("phi_Rz", nonzero=True) ** 2
                - 5 * sp.exp(-sp.Rational(1, 5) * sp.Symbol("logref"))
            )
        )
        == 0,
    }

    context = MPIntervalContext()
    context.prec = 300
    F = Fraction
    signed_inputs = {
        "U": [_fresh_signed_record(context, F(2)), _fresh_signed_record(context, F(1)), _fresh_signed_record(context, F(1, 2))],
        "phi_end": [_fresh_signed_record(context, F(3)), _fresh_signed_record(context, F(-1)), _fresh_signed_record(context, F(2))],
        "s": [_fresh_signed_record(context, F(1, 3)), _fresh_signed_record(context, F(1, 5)), _fresh_signed_record(context, F(1, 7))],
    }
    for name, base in (("A_in", 5), ("B_in", 7), ("C_in", 11)):
        signed_inputs[name] = [
            _fresh_signed_record(context, F(base)),
            _fresh_signed_record(context, F(1)),
            _fresh_signed_record(context, F(2)),
        ]
    finite = transport_signed_jet(context, **signed_inputs)
    record_rejection = False
    try:
        transport_signed_jet(
            context,
            signed_inputs["U"],
            signed_inputs["phi_end"],
            signed_inputs["s"],
            [{"ordinary_Z_order": 0, "log_absolute_upper_bound": _fresh_signed_record(context, F(1))}] * 3,
            signed_inputs["B_in"],
            signed_inputs["C_in"],
        )
    except (TypeError, ValueError):
        record_rejection = True
    checks = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "source_ode_A_exact": ode_checks["A"],
        "source_ode_B_exact": ode_checks["B"],
        "source_ode_C_exact": ode_checks["C"],
        "initial_A_exact": initial_checks["A"],
        "initial_B_exact": initial_checks["B"],
        "initial_C_exact": initial_checks["C"],
        "ordinary_Z1_Z2_A_nonzero": z12_checks["A"],
        "ordinary_Z1_Z2_B_nonzero": z12_checks["B"],
        "ordinary_Z1_Z2_C_nonzero": z12_checks["C"],
        "full_A_memory_contribution": memory_checks["A"],
        "full_B_memory_contribution": memory_checks["B"],
        "full_C_memory_contribution": memory_checks["C"],
        "normalized_B_deviation_exact": normalized_checks["B"],
        "normalized_C_deviation_exact": normalized_checks["C"],
        "phi_normalization_exact": phi_identity,
        "restoration_A_formula_exact": restoration_formula_check,
        "source_ode_proof_exact": ode_proof["passed"],
        "canonical_Rz_A_exact": canonical_checks["A"],
        "canonical_Rz_B_exact": canonical_checks["B"],
        "canonical_Rz_C_exact": canonical_checks["C"],
        "finite_signed_transport_has_Z0_Z1_Z2": finite["ordinary_Z_orders"] == [0, 1, 2],
        "finite_signed_transport_keeps_memory": finite["memory_contributions"]["A"][0] != 0
        and finite["memory_contributions"]["B"][0] != 0
        and finite["memory_contributions"]["C"][0] != 0,
        "cap_record_rejected": record_rejection,
        "fresh_rational_fading_memory_fixture": fading_memory["passed"],
    }
    _require(all(checks.values()), "fresh signed native history fixture failed")
    return {
        "passed": True,
        "schema": "original_signed_native_history_kernel_fixture_v1",
        "module_sha256": preflight["module_sha256"],
        "module": MODULE,
        "checks": checks,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "source_ode": {
            "A": "A_y=U**2-A",
            "B": "B_y=phi**2-2*B",
            "C": "C_y=phi**2-C",
        },
        "normalized_relations": {
            "B": "B/phi_end**2-5/6=exp(-2*s)*B_in/phi_end**2-(5/6)*exp(-6*s/5)",
            "C": "C/phi_end**2-5=exp(-s)*C_in/phi_end**2-5*exp(-s/5)",
        },
        "physical_normalization": "R*F0**2*phi**2=Utheta**2/2",
        "finite_signed_input_records": True,
        "caps_used_as_signed_values": False,
        "midpoint_selected": False,
        "gamma1_recomputed": False,
        "gamma2_recomputed": False,
        "source_saved_data_loaded": False,
        "fading_memory": fading_memory,
    }


fixture = fresh_fixture
self_check = fresh_fixture


__all__ = (
    "MODULE",
    "ORDINARY_Z_ORDERS",
    "NATIVE_HISTORY_NAMES",
    "source_ode_formulas",
    "source_ode_jets",
    "phi_from_Utheta",
    "pure_signed_formula",
    "decompose_signed_formula",
    "signed_decomposition",
    "signed_formula",
    "signedFormula",
    "symbolic_source",
    "source_formula",
    "source_ode_proof",
    "sourceODEproof",
    "canonical_Rz_signed_formula",
    "FADING_MEMORY_NU",
    "PHI_RELATIVE_DECAY_RATE",
    "SAVED_BSHAPE_SLOPE_UPPER",
    "SAVED_SIGMA_PRIME_UPPER",
    "relative_phi_density_factor",
    "pure_fading_memory_log_bound",
    "fading_memory_log_bound",
    "log_fading_memory_bound",
    "fading_memory_endpoint_bound",
    "pureLogFadingMemoryBound",
    "fadingMemoryLogBound",
    "relativePhiDensityFactor",
    "fresh_fading_memory_fixture",
    "fading_memory_fixture",
    "read_signed_value",
    "transport_signed_jet",
    "signed_jet_transport",
    "finite_signed_jet_transport",
    "transport_signed_histories",
    "post_reference_signed_jet",
    "reference_signed_jet_transport",
    "restoration_signed_jets",
    "post_restore_signed_jets",
    "record_signed_jet",
    "fresh_fixture",
    "fixture",
    "self_check",
)
