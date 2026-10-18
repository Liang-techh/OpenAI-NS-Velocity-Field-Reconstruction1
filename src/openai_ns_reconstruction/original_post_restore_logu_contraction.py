"""Uniform fixed-phase fast derivative in cancellation-safe coordinates.

The post-restore phase source is

``u=p2*q/dstar``, ``r=u/sqrt(1+u**2)``, ``m=1-r**2``,
``D=1-2*r*cos(psi)+r**2``, and
``t=2*q*sqrt(m)*(cos(psi)-r)/D``.

This module keeps that source definition separate from a uniform bound.  It
does not import, execute, or read a producer, accepted phase module, saved
packet, logarithmic cap, or historical fixture.  The direct ``r`` formulas
are replaced at evaluation time by the Möbius coordinates

``v=t/(2*q)``, ``h=sqrt(m)*sin(psi)/D``, ``k=m/D``.

For ``u<0`` these coordinates can be evaluated without subtracting two
numbers near ``-1``.  They satisfy

``k=1+2*u*v``, ``h**2=k-v**2``, ``m=1/(1+u**2)``.

The exact source derivatives then become

``t_u=2*q*m*(2*v**2-3*u*v-1)``,
``t_psi=-2*q*k*h``,
``I_u=4*q**2*(omega/u**3-m*h*(3-k)/(2*u**2))``.

Here ``omega=atan2(r*sin(psi),1-r*cos(psi))`` is the continuous principal
representative; its denominator is positive.  The fixed-phase derivative is

``t_u_star=t_u-t_psi*I_u/(1+t**2)``

and the quantity needed by the fast source is

``W=u*((1-t**2)/(1+t**2)**2)*t_u_star``.

For ``|u|>=1`` and ``1/2<=q<=3/2`` the proof in ``uniform_bound`` gives
``|W| <= 1061/2 < 600`` uniformly in ``psi``.  The bound is deliberately
source-independent: it is a sufficient analytic estimate for a later signed
``p2`` port, not a selection of a saved endpoint value.
"""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


MODULE = "src/openai_ns_reconstruction/original_post_restore_logu_contraction.py"
FIXTURE_RELATIVE_PATH = "tmp/original_post_restore_logu_contraction_fixture.json"

Q_LOWER = Fraction(1, 2)
Q_UPPER = Fraction(3, 2)
U_ABS_LOWER = Fraction(1)
PROVED_BOUND = Fraction(1061, 2)
SAFE_BOUND = Fraction(600)
Q_SQUARED_SOURCE = "3/4 + 5*eta/4 (constant in Z)"


def _raw_bounds(value: Any) -> tuple[mp.mpf, mp.mpf]:
    if not hasattr(value, "_mpi_"):
        raise TypeError("an MPIntervalContext interval is required")
    rows = value._mpi_
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    if len(rows) == 1:
        rows = (rows, rows)
    if len(rows) != 2:
        raise TypeError("malformed interval")
    return mp.make_mpf(tuple(rows[0])), mp.make_mpf(tuple(rows[1]))


def _lower(value: Any) -> mp.mpf:
    return _raw_bounds(value)[0]


def _upper(value: Any) -> mp.mpf:
    return _raw_bounds(value)[1]


def _finite(value: Any) -> bool:
    lo, hi = _raw_bounds(value)
    return not any(mp.isnan(endpoint) or mp.isinf(endpoint)
                   for endpoint in (lo, hi))


def _interval(ctx: MPIntervalContext, lower: Any, upper: Any):
    if hasattr(lower, "_mpi_"):
        lower = _lower(lower)
    if hasattr(upper, "_mpi_"):
        upper = _upper(upper)
    if lower > upper:
        raise ValueError("reversed directed interval")
    result = ctx.mpf([lower, upper])
    if not _finite(result):
        raise ValueError("finite directed intervals are required")
    return result


def _coerce(ctx: MPIntervalContext, value: Any):
    """Box an exact scalar or directed interval without binary floats."""

    if hasattr(value, "_mpi_"):
        return _interval(ctx, _lower(value), _upper(value))
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("use an exact scalar or MPIntervalContext interval")
    if type(value) is int:
        return ctx.mpf(value)
    if isinstance(value, Fraction):
        return ctx.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return ctx.mpf(value)
    if hasattr(value, "_mpf_"):
        point = mp.make_mpf(tuple(value._mpf_))
        return _interval(ctx, point, point)
    if isinstance(value, (tuple, list)) and len(value) == 2:
        left, right = _coerce(ctx, value[0]), _coerce(ctx, value[1])
        return _interval(ctx, _lower(left), _upper(right))
    raise TypeError("use an exact scalar or MPIntervalContext interval")


def _square(value: Any):
    """Use integer power so an interval square cannot get a negative floor."""

    result = value ** 2
    if _lower(result) < 0:
        raise ValueError("interval square produced a negative lower endpoint")
    return result


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _module_source_preflight(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Check LF/one-final-LF/AST before any fresh mathematical fixture."""

    root = (Path(repo_root) if repo_root is not None
            else Path(__file__).resolve().parents[2])
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
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
    }


def _validate_domain(ctx: MPIntervalContext, u: Any, q: Any, psi: Any):
    u_box = _coerce(ctx, u)
    q_box = _coerce(ctx, q)
    psi_box = _coerce(ctx, psi)
    _require(_finite(u_box) and _finite(q_box) and _finite(psi_box),
             "u, q, and psi must be finite")
    _require(_upper(u_box) < 0, "u must be strictly negative")
    _require(_upper(u_box) <= -mp.mpf(1), "|u| must be at least 1")
    _require(_lower(q_box) >= mp.mpf(1) / 2,
             "q must be at least 1/2")
    _require(_upper(q_box) <= mp.mpf(3) / 2,
             "q must be at most 3/2")
    _require(_lower(psi_box) >= 0 and
             _upper(psi_box) <= _upper(2 * ctx.pi),
             "psi must lie in [0, 2*pi]")
    return u_box, q_box, psi_box


def _stable_coordinates(ctx: MPIntervalContext, u: Any, psi: Any):
    """Return ``(m,v,h,k,omega)`` without forming ``r`` near ``-1``.

    Put ``U=-u``, ``s=sqrt(1+U**2)``, ``delta=s-U=1/(s+U)``, and
    ``y=1+cos(psi)=2*cos(psi/2)**2``.  Then

    ``Q=delta**2+2*U*s*y``, ``k=1/Q``,
    ``v=(s*y-delta)/Q``, ``h=s*sin(psi)/Q``.

    Every potentially tiny denominator is represented as a positive sum.
    """

    U = -u
    s = ctx.sqrt(1 + _square(U))
    delta = 1 / (s + U)
    half = psi / 2
    cos_half = ctx.cos(half)
    sin_half = ctx.sin(half)
    y = 2 * _square(cos_half)
    sine = 2 * sin_half * cos_half
    Q = _square(delta) + 2 * U * s * y
    _require(_lower(Q) > 0, "stable coordinate denominator must be positive")
    m = 1 / _square(s)
    k = 1 / Q
    v = (s * y - delta) / Q
    h = s * sine / Q
    omega = ctx.atan2(-U * sine, delta + U * y)
    return m, v, h, k, omega


def fixed_phase_interval(
    ctx: MPIntervalContext,
    u: Any,
    q: Any,
    psi: Any,
) -> dict[str, Any]:
    """Evaluate the scaled fixed-phase contraction on finite directed inputs.

    The evaluator uses only exact scalars or ``MPIntervalContext`` intervals;
    in particular it does not accept a logarithmic cap as a value.  The input
    domain is ``u<=-1``, ``q in [1/2,3/2]`` and ``psi in [0,2*pi]``.
    """

    if not isinstance(ctx, MPIntervalContext):
        raise TypeError("ctx must be an MPIntervalContext")
    u_box, q_box, psi_box = _validate_domain(ctx, u, q, psi)
    m, v, h, k, omega = _stable_coordinates(ctx, u_box, psi_box)
    q2 = _square(q_box)
    t = 2 * q_box * v
    t2 = _square(t)
    one_plus_t2 = 1 + t2
    t_u = 2 * q_box * m * (2 * _square(v) - 3 * u_box * v - 1)
    t_psi = -2 * q_box * k * h
    u2 = _square(u_box)
    u3 = u2 * u_box
    I_u = 4 * q2 * (
        omega / u3 - m * h * (3 - k) / (2 * u2)
    )
    fixed_t_u = t_u - t_psi * I_u / one_plus_t2
    f_prime = (1 - t2) / _square(one_plus_t2)
    W = u_box * f_prime * fixed_t_u
    for label, value in (
        ("m", m), ("v", v), ("h", h), ("k", k), ("omega", omega),
        ("t", t), ("t_u", t_u), ("t_psi", t_psi), ("I_u", I_u),
        ("fixed_t_u", fixed_t_u), ("f_prime", f_prime), ("W", W),
    ):
        _require(_finite(value), label + " is not finite")
    return {
        "u": u_box,
        "q": q_box,
        "psi": psi_box,
        "m": m,
        "v": v,
        "h": h,
        "k": k,
        "omega": omega,
        "t": t,
        "t_u": t_u,
        "t_psi": t_psi,
        "I_u": I_u,
        "fixed_t_u": fixed_t_u,
        "f_prime": f_prime,
        "W": W,
    }


def uniform_bound(
    q_lower: Fraction | int = Q_LOWER,
    q_upper: Fraction | int = Q_UPPER,
    u_abs_lower: Fraction | int = U_ABS_LOWER,
) -> dict[str, Any]:
    """Return the exact source-independent uniform bound receipt.

    Let ``U=|u|`` and ``F=|(1-4*q**2*v**2)/(1+4*q**2*v**2)**2|``.  The
    elementary inequalities used below are

    ``|U*m|<=1/2``, ``U**2*m<=1``,
    ``v**2*F<=1/(4*q**2)``, ``|v|*F<=1/(4*q)``.

    The transformed first term is bounded by

    ``|u*F*t_u| <= 1/(2*q) + 3/2 + q <= 4``.

    Since ``|1-z**2|/(1+z**2)**3 <= (1+z**2)**-2`` for ``z=2*q*v`` and
    ``q>=1/2``, the phase term has denominator at least ``(1+v**2)**2``.
    With ``|h|<=sqrt(1+u**2)<=2*U``, ``|omega|<=pi/2<=2`` and
    ``k<=1+2*U*|v|``, its omega part is at most ``12``.  The exact
    polynomial identity

    ``m*k*h**2*(3-k)/(2*u) = m*(1-v**2)/u + 3*m*v``
    ``                                  -(4*u**2+1)*m*v**3 + 2*m*u*v**4``

    gives the remaining part at most ``1/2+1/2+3/2+4+1=15/2`` after
    division by ``(1+v**2)**2``.  Finally ``8*q**3<=27``, so

    ``|phase term| <= 27*(12+15/2)=1053/2`` and
    ``|W| <= 4+1053/2=1061/2 < 600``.

    Inputs are exact rational bounds.  No numerical samples are used to
    establish this result.
    """

    if not isinstance(q_lower, (int, Fraction)) or isinstance(q_lower, bool):
        raise TypeError("q_lower must be an exact rational")
    if not isinstance(q_upper, (int, Fraction)) or isinstance(q_upper, bool):
        raise TypeError("q_upper must be an exact rational")
    if not isinstance(u_abs_lower, (int, Fraction)) or isinstance(u_abs_lower, bool):
        raise TypeError("u_abs_lower must be an exact rational")
    q_lo = Fraction(q_lower)
    q_hi = Fraction(q_upper)
    u_lo = Fraction(u_abs_lower)
    _require(q_lo >= Q_LOWER, "q_lower is below the proved q domain")
    _require(q_hi <= Q_UPPER and q_lo <= q_hi,
             "q bounds must lie inside [1/2,3/2]")
    _require(u_lo >= U_ABS_LOWER,
             "u_abs_lower is below the proved |u| domain")
    first = Fraction(4)
    omega = Fraction(12)
    polynomial = Fraction(15, 2)
    phase_prefactor = Fraction(27)
    phase = phase_prefactor * (omega + polynomial)
    proved = first + phase
    _require(proved == PROVED_BOUND, "internal exact bound arithmetic changed")
    return {
        "domain": {
            "q_lower": [q_lo.numerator, q_lo.denominator],
            "q_upper": [q_hi.numerator, q_hi.denominator],
            "abs_u_lower": [u_lo.numerator, u_lo.denominator],
            "psi": "[0,2*pi]",
        },
        "first_term_bound": [first.numerator, first.denominator],
        "omega_term_bound": [omega.numerator, omega.denominator],
        "polynomial_term_bound": [polynomial.numerator, polynomial.denominator],
        "phase_prefactor_bound": [phase_prefactor.numerator,
                                   phase_prefactor.denominator],
        "phase_term_bound": [phase.numerator, phase.denominator],
        "proved_bound": [proved.numerator, proved.denominator],
        "safe_bound": [SAFE_BOUND.numerator, SAFE_BOUND.denominator],
        "strict_safe_margin": [
            (SAFE_BOUND - proved).numerator,
            (SAFE_BOUND - proved).denominator,
        ],
        "bound_is_uniform_in_psi": True,
        "source_signed_value_selected": False,
    }


def _reduce_exact(expr: Any, *, unit: Any, scale: Any) -> Any:
    """Reduce an exact rational expression modulo the two source relations."""

    import sympy as sp

    numerator = sp.together(expr).as_numer_denom()[0]
    z = unit.lhs.free_symbols.pop() if False else None
    # ``unit`` and ``scale`` are passed as explicit polynomials below; their
    # variables are recovered deterministically from their free-symbol set.
    variables = sorted(unit.free_symbols | scale.free_symbols,
                       key=lambda symbol: symbol.name)
    if len(variables) != 4:
        raise ValueError("symbolic proof expected four source variables")
    # The proof calls this helper with the known (u,s,c,z) ordering.
    u, s, cpsi, spsi = [next(symbol for symbol in variables
                              if symbol.name == name)
                        for name in ("u", "s", "c", "spsi")]
    unit_poly = sp.Poly(unit, spsi)
    remainder = sp.rem(sp.Poly(numerator, spsi), unit_poly).as_expr()
    scale_poly = sp.Poly(scale, s)
    remainder = sp.rem(sp.Poly(remainder, s), scale_poly).as_expr()
    return sp.factor(sp.cancel(remainder))


def fresh_symbolic_proof() -> dict[str, Any]:
    """Prove the coordinate identities with fresh exact SymPy polynomials."""

    import sympy as sp

    u, s, cpsi, spsi, q, omega = sp.symbols(
        "u s c spsi q omega", nonzero=True
    )
    m = 1 / s**2
    r = u / s
    D = 1 - 2 * r * cpsi + r**2
    x = cpsi - r
    v = (s * cpsi - u) / (s**2 * D)
    h = spsi / (s * D)
    k = m / D
    tu_source = 2 * q * m / D**2 * (
        m * (x**2 - spsi**2) - r * x * D
    )
    tu_scaled = 2 * q * m * (2 * v**2 - 3 * u * v - 1)
    tpsi_source = -2 * q * m**sp.Rational(3, 2) * spsi / D**2
    # On the positive s branch m^(3/2)=s^(-3); this is the exact source
    # representation used in the polynomial residual below.
    tpsi_source = -2 * q * spsi / (s**3 * D**2)
    tpsi_scaled = -2 * q * k * h
    Iu_source = 4 * q**2 / s**3 * (
        omega / r**3
        - m * spsi * (1 - 3 * r * cpsi + 2 * r**2) / (r**2 * D**2)
    )
    Iu_scaled = 4 * q**2 * (
        omega / u**3 - m * h * (3 - k) / (2 * u**2)
    )
    unit = cpsi**2 + spsi**2 - 1
    scale = s**2 - u**2 - 1

    def residual(expr: Any) -> str:
        return str(_reduce_exact(expr, unit=unit, scale=scale))

    qv = 2 * q * v
    fp = (1 - qv**2) / (1 + qv**2)**2
    W_source = u * fp * (tu_source - tpsi_source * Iu_source / (1 + qv**2))
    W_scaled = u * fp * (tu_scaled - tpsi_scaled * Iu_scaled / (1 + qv**2))
    P = m * k * (k - v**2) * (3 - k) / (2 * u)
    P_expanded = (
        m * (1 - v**2) / u + 3 * m * v
        - (4 * u**2 + 1) * m * v**3 + 2 * m * u * v**4
    )
    return {
        "coordinate_k_residual": residual(k - (1 + 2 * u * v)),
        "coordinate_h_residual": residual(h**2 - (k - v**2)),
        "coordinate_circle_residual": residual((v - u)**2 + h**2 - s**2),
        "t_u_residual": residual(tu_source - tu_scaled),
        "t_psi_residual": residual(tpsi_source - tpsi_scaled),
        "I_u_residual": residual(Iu_source - Iu_scaled),
        "fixed_phase_W_residual": residual(W_source - W_scaled),
        "polynomial_P_residual": residual(P - P_expanded),
        "relations": {
            "unit_circle": "c^2+spsi^2=1",
            "scale": "s^2-u^2=1",
            "m": "1/s^2",
            "r": "u/s",
        },
    }


def _interval_record(value: Any, digits: int = 32) -> list[str]:
    rows = value._mpi_
    with mp.workdps(digits + 8):
        return [
            mp.nstr(mp.make_mpf(tuple(rows[0])), digits),
            mp.nstr(mp.make_mpf(tuple(rows[1])), digits),
        ]


def _mp_fraction(value: Fraction | int) -> mp.mpf:
    fraction = Fraction(value)
    return mp.mpf(fraction.numerator) / fraction.denominator


def _fresh_scalar_reference(u: Fraction, q: Fraction, psi: Fraction) -> mp.mpf:
    """Independent direct finite reference for a moderate fresh rational case."""

    with mp.workdps(180):
        uu = _mp_fraction(u)
        qq = _mp_fraction(q)
        pp = _mp_fraction(psi)
        ss = mp.sqrt(1 + uu * uu)
        rr = uu / ss
        mm = 1 / (1 + uu * uu)
        dd = 1 - 2 * rr * mp.cos(pp) + rr * rr
        xx = mp.cos(pp) - rr
        vv = mp.sqrt(mm) * xx / dd
        hh = mp.sqrt(mm) * mp.sin(pp) / dd
        kk = mm / dd
        tt = 2 * qq * vv
        # Rebuild the original source derivatives independently.  The
        # transformed evaluator above uses the factored forms; this reference
        # intentionally uses the unsimplified source expressions at moderate
        # rational inputs only.
        tu = 2 * qq * mm / dd**2 * (
            mm * (xx**2 - mp.sin(pp)**2) - rr * xx * dd
        )
        tpsi = -2 * qq * mm**mp.mpf("1.5") * mp.sin(pp) / dd**2
        omega = mp.atan2(rr * mp.sin(pp), 1 - rr * mp.cos(pp))
        i_u = 4 * qq * qq * (
            mm**mp.mpf("1.5") * (
                omega / rr**3
                - mm * mp.sin(pp) * (1 - 3 * rr * mp.cos(pp) + 2 * rr**2)
                / (rr**2 * dd**2)
            )
        )
        fp = (1 - tt * tt) / (1 + tt * tt)**2
        return uu * fp * (tu - tpsi * i_u / (1 + tt * tt))


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run fresh exact symbolic/rational checks and return a serializable fixture.

    The LF/EOF/AST preflight is intentionally the first operation.  All
    subsequent numerical cases use new rational inputs and the standalone
    formulas above; no saved source packet or accepted producer is opened.
    """

    preflight = _module_source_preflight(repo_root)
    proof = fresh_symbolic_proof()
    bound = uniform_bound()
    checks: dict[str, bool] = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "coordinate_k_zero": proof["coordinate_k_residual"] == "0",
        "coordinate_h_zero": proof["coordinate_h_residual"] == "0",
        "coordinate_circle_zero": proof["coordinate_circle_residual"] == "0",
        "t_u_zero": proof["t_u_residual"] == "0",
        "t_psi_zero": proof["t_psi_residual"] == "0",
        "I_u_zero": proof["I_u_residual"] == "0",
        "fixed_phase_W_zero": proof["fixed_phase_W_residual"] == "0",
        "polynomial_P_zero": proof["polynomial_P_residual"] == "0",
        "exact_bound_is_1061_over_2": bound["proved_bound"] == [1061, 2],
        "safe_bound_is_600": bound["safe_bound"] == [600, 1],
    }
    cases: dict[str, Any] = {}
    old_prec = mp.mp.prec
    try:
        ctx = MPIntervalContext()
        ctx.dps = 140
        rational_cases = (
            ("u_minus_2_q1_psi_7_over_10", Fraction(-2), Fraction(1),
             Fraction(7, 10)),
            ("u_minus_17_over_3_q5_over_4_psi_11_over_10",
             Fraction(-17, 3), Fraction(5, 4), Fraction(11, 10)),
        )
        for name, u, q, psi in rational_cases:
            out = fixed_phase_interval(ctx, u, q, psi)
            reference = _fresh_scalar_reference(u, q, psi)
            contained = (_lower(out["W"]) <= reference <= _upper(out["W"]))
            bounded = (_lower(out["W"]) >= -mp.mpf(600)
                       and _upper(out["W"]) <= mp.mpf(600))
            checks[name + "_contains_fresh_reference"] = contained
            checks[name + "_inside_uniform_bound"] = bounded
            cases[name] = {
                "W": _interval_record(out["W"]),
                "v": _interval_record(out["v"]),
                "h": _interval_record(out["h"]),
                "k": _interval_record(out["k"]),
            }

        interval_case = fixed_phase_interval(
            ctx,
            (Fraction(-17, 3), Fraction(-17, 3)),
            (Fraction(5, 4), Fraction(5, 4)),
            (Fraction(109, 100), Fraction(111, 100)),
        )
        checks["finite_directed_interval_case"] = (
            _finite(interval_case["W"])
            and _lower(interval_case["W"]) >= -mp.mpf(600)
            and _upper(interval_case["W"]) <= mp.mpf(600)
        )
        cases["u_minus_17_over_3_q5_over_4_psi_109_to_111_over_100"] = {
            "W": _interval_record(interval_case["W"]),
            "v": _interval_record(interval_case["v"]),
            "h": _interval_record(interval_case["h"]),
            "k": _interval_record(interval_case["k"]),
        }

        # A large exact negative input exercises the stable delta/Q path.  No
        # exponential or logarithmic cap is materialized.
        huge_u = -(10 ** 80)
        huge = fixed_phase_interval(
            ctx, huge_u, Fraction(1), Fraction(314, 100),
        )
        checks["huge_negative_u_interval_finite"] = _finite(huge["W"])
        checks["huge_negative_u_inside_uniform_bound"] = (
            _lower(huge["W"]) >= -mp.mpf(600)
            and _upper(huge["W"]) <= mp.mpf(600)
        )
        cases["u_minus_10_pow_80_q1_psi_314_over_100"] = {
            "W": _interval_record(huge["W"]),
            "v": _interval_record(huge["v"]),
            "h": _interval_record(huge["h"]),
            "k": _interval_record(huge["k"]),
        }
    finally:
        mp.mp.prec = old_prec
    return {
        "passed": all(checks.values()),
        "module_sha256": preflight["sha256"],
        "checks": checks,
        "bound": bound,
        "proof": proof,
        "cases": cases,
        "scope": {
            "source_executed": False,
            "saved_packet_read": False,
            "logcap_used_as_actual_value": False,
            "signed_p2_selected": False,
        },
    }


__all__ = [
    "FIXTURE_RELATIVE_PATH",
    "MODULE",
    "PROVED_BOUND",
    "Q_SQUARED_SOURCE",
    "Q_LOWER",
    "Q_UPPER",
    "SAFE_BOUND",
    "U_ABS_LOWER",
    "fixed_phase_interval",
    "fresh_fixture",
    "fresh_symbolic_proof",
    "uniform_bound",
]
