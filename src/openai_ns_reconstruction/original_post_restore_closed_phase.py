"""Pure directed closed phase source for the post-restore fast loop.

The source is the post-restore specialization of the checked C5 loop:

``t = 2*q/s*(cos(psi)-r)/(1-2*r*cos(psi)+r**2)``

with ``s=sqrt(1+u**2)`` and ``r=u/s``.  This module evaluates only the
closed primitive and the fixed-phase signed contraction.  It does not import,
execute, or read any producer, owner, checker, recurrence, quadrature, saved
packet, or logarithmic cap.

Inputs are exact scalars or finite directed ``MPIntervalContext`` intervals.
The formulas use ``atan2`` because the interval context exposes ``atan2`` but
not ``atan``.  The denominator of that ``atan2`` is checked to be strictly
positive, so the returned angle is the continuous principal representative
for ``0 <= psi <= 2*pi``.

The ``r=0`` expression is handled by its exact removable limit.  A
non-singleton interval containing zero is rejected before any ``1/r**2`` or
``1/r**3`` expression is formed.  The explicit ``psi="2*pi"`` token selects
the named upper endpoint identity; an ordinary interval such as ``2*c.pi``
keeps interval trigonometry and is widened only enough to contain that named
endpoint.  The module binds its formula to the read-only C5 ``loop_functions``
source by both byte SHA256 and function AST SHA256; this is source provenance
only and does not bind actual signed source parameters.
"""

from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


MODULE = "src/openai_ns_reconstruction/original_post_restore_closed_phase.py"
FIXTURE_RELATIVE_PATH = "tmp/post_restore_closed_phase_fixture.json"

C5_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
C5_SOURCE_SHA256 = (
    "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
)
C5_LOOP_FUNCTION = "loop_functions"
C5_LOOP_FUNCTION_AST_SHA256 = (
    "fffc175d358e7cce71d3046ec10303e22d46b5ebb3145f9a44c0d1fe88796f8c"
)
C5_LOOP_OUTPUT_PORTS = (
    "B_over_Pstar_phi",
    "direction_at_inverse",
    "phase_lambda_at_inverse",
)

# This module is formula-only.  No actual signed ``u`` or ``q`` source value
# is selected from a saved report.
SIGNED_SOURCE_ACTUAL_PARAMETERS_BOUND = False


def _raw_bounds(value: Any) -> tuple[mp.mpf, mp.mpf]:
    if not hasattr(value, "_mpi_"):
        raise TypeError("An MPIntervalContext interval is required")
    rows = value._mpi_
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    if len(rows) == 1:
        rows = (rows, rows)
    if len(rows) != 2:
        raise TypeError("Malformed interval")
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
        raise ValueError("Reversed directed interval")
    result = ctx.mpf([lower, upper])
    if not _finite(result):
        raise ValueError("finite directed intervals are required")
    return result


def _coerce(ctx: MPIntervalContext, value: Any):
    """Box an exact scalar or directed interval without binary floats."""

    if hasattr(value, "_mpi_"):
        return _interval(ctx, _lower(value), _upper(value))
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("Use an exact scalar or MPIntervalContext interval")
    if type(value) is int:
        return ctx.mpf(value)
    if isinstance(value, Fraction):
        return ctx.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return ctx.mpf(value)
    if hasattr(value, "_mpf_"):
        # ``mp.mpf(value)`` obeys the process-global precision and can round
        # a high-precision scalar to 53 bits before the IV context sees it.
        # Rebuild directly from the exact MPF tuple instead.
        point = mp.make_mpf(tuple(value._mpf_))
        return _interval(ctx, point, point)
    if isinstance(value, (tuple, list)):
        if len(value) == 4 and all(type(part) is int for part in value):
            point = mp.make_mpf(tuple(value))
            return _interval(ctx, point, point)
        if len(value) == 2:
            left, right = _coerce(ctx, value[0]), _coerce(ctx, value[1])
            return _interval(ctx, _lower(left), _upper(right))
    raise TypeError("Use an exact scalar or MPIntervalContext interval")


def _singleton(value: Any) -> bool:
    return _lower(value) == _upper(value)


def _square(ctx: MPIntervalContext, value: Any):
    """Safe square, including the zero-crossing case.

    ``MPIntervalContext``'s binary multiplication can retain a negative lower
    endpoint for ``x*x`` when ``x`` crosses zero.  Integer power has the
    intended interval semantics, so use it for every square in this module.
    """

    result = value ** 2
    if _lower(result) < 0:
        raise ValueError("interval square produced a negative lower endpoint")
    return result


def _require_positive(value: Any, label: str) -> None:
    if _lower(value) <= 0:
        raise ValueError(label + " must be strictly positive on the input interval")


def _hull(ctx: MPIntervalContext, *values: Any):
    if not values:
        raise ValueError("an interval hull needs at least one value")
    return _interval(
        ctx,
        min(_lower(value) for value in values),
        max(_upper(value) for value in values),
    )


def _require_angle(ctx: MPIntervalContext, psi: Any) -> None:
    if _lower(psi) < 0:
        raise ValueError("psi must lie in [0, 2*pi]")
    two_pi = 2 * ctx.pi
    if _upper(psi) > _upper(two_pi):
        raise ValueError("psi must lie in [0, 2*pi]")


def c5_loop_source_binding(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Return and verify the read-only C5 source byte/AST binding."""

    root = (Path(repo_root) if repo_root is not None
            else Path(__file__).resolve().parents[2])
    path = root / Path(C5_SOURCE)
    raw = path.read_bytes()
    tree = ast.parse(raw.decode("utf-8"), filename=str(path))
    functions = [
        node for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == C5_LOOP_FUNCTION
    ]
    if len(functions) != 1:
        raise ValueError("expected exactly one original C5 loop_functions AST")
    function = functions[0]
    source_sha = hashlib.sha256(raw).hexdigest()
    ast_sha = hashlib.sha256(
        ast.dump(function, include_attributes=False).encode("utf-8")
    ).hexdigest()
    if source_sha != C5_SOURCE_SHA256 or ast_sha != C5_LOOP_FUNCTION_AST_SHA256:
        raise ValueError("original C5 loop_functions source binding changed")
    return {
        "path": C5_SOURCE,
        "sha256": source_sha,
        "function": C5_LOOP_FUNCTION,
        "function_ast_sha256": ast_sha,
        "function_line_span": [function.lineno, function.end_lineno],
        "output_ports": list(C5_LOOP_OUTPUT_PORTS),
        "source_read_only": True,
        "source_executed": False,
        "dag_reconstructed": False,
        "signed_source_actual_parameters_bound": SIGNED_SOURCE_ACTUAL_PARAMETERS_BOUND,
    }


def closed_phase(
    c: MPIntervalContext,
    u: Any,
    q: Any,
    psi: Any,
    *,
    a: Any = Fraction(4, 5),
) -> dict[str, Any]:
    """Evaluate the post-restore closed phase source on directed intervals.

    ``u``, ``q`` and ``psi`` are finite signed intervals.  The phase angle is
    restricted to ``[0, 2*pi]``.  A non-singleton ``r`` interval containing
    zero is rejected because the direct primitive uses removable ``r**-2``
    and ``r**-3`` terms; the exact singleton ``r=0`` branch uses its closed
    limits instead.

    The returned ``contractedC`` is

    ``(1+t**2)*t_u - t_psi*I_u``.

    ``bphi_u_over_e`` is the normalized signed source derivative requested by
    the post-restore row, with the exact API parameter ``a`` defaulting to
    ``4/5``.  No actual source value is inferred from a log cap.
    """

    if not isinstance(c, MPIntervalContext):
        raise TypeError("c must be an MPIntervalContext")
    u_box = _coerce(c, u)
    q_box = _coerce(c, q)
    # A string token is the explicit typed endpoint form.  Ordinary IV
    # expressions, including ``2*c.pi``, remain ordinary intervals.
    psi_token = psi.strip().lower() if isinstance(psi, str) else None
    if psi_token == "0":
        psi_named_endpoint = "zero"
        psi_box = c.mpf(0)
    elif psi_token == "2*pi":
        psi_named_endpoint = "two_pi"
        psi_box = 2 * c.pi
    else:
        psi_named_endpoint = None
        psi_box = _coerce(c, psi)
    a_box = _coerce(c, a)
    for label, value in (("u", u_box), ("q", q_box),
                         ("psi", psi_box), ("a", a_box)):
        if not _finite(value):
            raise ValueError(label + " must be finite")
    _require_angle(c, psi_box)

    u2 = _square(c, u_box)
    q2 = _square(c, q_box)
    one = c.mpf(1)
    s = c.sqrt(one + u2)
    m = one / (one + u2)
    r = u_box / s
    r_zero = _lower(r) == 0 and _upper(r) == 0
    if _lower(r) <= 0 <= _upper(r) and not r_zero:
        raise ValueError(
            "r interval contains zero without being singleton; "
            "use a narrower directed interval"
        )

    two_pi = 2 * c.pi
    at_zero = psi_named_endpoint == "zero" or (
        psi_named_endpoint is None
        and _lower(psi_box) == 0
        and _upper(psi_box) == 0
    )
    at_two_pi = psi_named_endpoint == "two_pi"
    named_endpoint = at_zero or at_two_pi
    ordinary_upper_endpoint = (
        not named_endpoint
        and _lower(psi_box) <= _lower(two_pi)
        and _upper(psi_box) >= _upper(two_pi)
    )
    if named_endpoint:
        # Only the explicit token (or exact singleton zero) selects these
        # endpoint identities.  A regular interval never gets canonicalized.
        cpsi = c.mpf(1)
        spsi = c.mpf(0)
        omega = c.mpf(0)
    else:
        cpsi = c.cos(psi_box)
        spsi = c.sin(psi_box)
        omega = c.atan2(r * spsi, one - r * cpsi)
    x = cpsi - r
    x2 = _square(c, x)
    y2 = _square(c, spsi)
    denominator = x2 + y2
    _require_positive(denominator, "D")
    one_minus_rc = one - r * cpsi
    _require_positive(one_minus_rc, "1-r*cos(psi)")

    alpha = 2 * q_box / s
    t = alpha * x / denominator
    t2 = _square(c, t)
    one_plus_t2 = one + t2
    _require_positive(one_plus_t2, "1+t^2")
    m32 = m * c.sqrt(m)

    if r_zero:
        # Exact removable r=0 limits.  These formulas remain interval-safe
        # when q or psi is a non-singleton interval, while u itself is exact 0.
        if named_endpoint:
            t_u = 2 * q_box
            t_psi = c.mpf(0)
            I = c.mpf(0) if at_zero else 4 * c.pi * q2
            I_u = c.mpf(0)
        else:
            t_u = 2 * q_box * c.cos(2 * psi_box)
            t_psi = -2 * q_box * spsi
            I = 2 * q2 * psi_box + q2 * c.sin(2 * psi_box)
            I_u = 4 * q2 * (
                spsi + c.sin(3 * psi_box) / 3
            )
        omega = c.mpf(0)
    else:
        r2 = _square(c, r)
        _require_positive(r2, "r^2")
        r3 = r2 * r
        D2 = _square(c, denominator)
        L = one - 3 * r * cpsi + 2 * r2
        t_u = (2 * q_box * m / D2) * (
            m * (x2 - y2) - r * x * denominator
        )
        if named_endpoint:
            t_psi = c.mpf(0)
            I = c.mpf(0) if at_zero else 4 * c.pi * q2
            I_u = c.mpf(0)
        else:
            t_psi = -2 * q_box * m32 * spsi / D2
            I = 2 * q2 * (
                psi_box
                + m * spsi / (r * denominator)
                + (3 * r2 - 1) * omega / r2
            )
            I_u = 4 * q2 * m32 * (
                omega / r3
                - m * spsi * L / (r2 * D2)
            )

    if ordinary_upper_endpoint:
        # The ordinary ``2*c.pi`` interval encloses the named upper endpoint,
        # but remains a genuine interval.  Hull only the known endpoint values
        # so cancellation in the primitive cannot exclude the exact identity.
        I = _hull(c, I, 4 * c.pi * q2)
        I_u = _hull(c, I_u, c.mpf(0))

    K = one / (two_pi * (one + 2 * q2))
    Phi = c.mpf(0) if at_zero else c.mpf(1) if at_two_pi else K * (psi_box + I)
    if ordinary_upper_endpoint:
        Phi = _hull(c, Phi, c.mpf(1))
    lam = K * one_plus_t2
    contracted = one_plus_t2 * t_u - t_psi * I_u
    fixedphase_t_u = contracted / one_plus_t2
    bphi_u_over_e = (
        -a_box * (one + 2 * q2) / 2
        * (one - t2) / (one_plus_t2 ** 3)
        * contracted
    )

    return {
        "t": t,
        "t_u": t_u,
        "t_psi": t_psi,
        "I": I,
        "I_u": I_u,
        "Phi": Phi,
        "lambda": lam,
        "contractedC": contracted,
        "fixedphase_t_u": fixedphase_t_u,
        "bphi_u_over_e": bphi_u_over_e,
        "r": r,
        "s": s,
        "m": m,
        "D": denominator,
        "omega": omega,
        "K": K,
        "a": a_box,
        "signed_source_actual_parameters_bound": SIGNED_SOURCE_ACTUAL_PARAMETERS_BOUND,
    }


def _module_source_preflight(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Check this module's LF/EOF/AST before any fresh numerical work."""

    root = (Path(repo_root) if repo_root is not None
            else Path(__file__).resolve().parents[2])
    path = root / Path(MODULE)
    raw = path.read_bytes()
    if b"\r" in raw:
        raise ValueError("module source must use LF line endings")
    if not raw.endswith(b"\n") or raw.endswith(b"\n\n"):
        raise ValueError("module source must have one final LF")
    ast.parse(raw.decode("utf-8"), filename=str(path))
    return {
        "passed": True,
        "path": MODULE,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
    }


def _interval_record(value: Any, digits: int = 140) -> list[str]:
    rows = value._mpi_
    with mp.workdps(digits + 10):
        return [
            mp.nstr(mp.make_mpf(tuple(rows[0])), digits),
            mp.nstr(mp.make_mpf(tuple(rows[1])), digits),
        ]


def _fresh_symbolic_proof() -> dict[str, Any]:
    """Return exact rational residuals for the nonzero-r primitive identities."""

    import sympy as sp

    r, cpsi, spsi, omega, psi, qv = sp.symbols(
        "r c s omega psi q", nonzero=True
    )
    D = 1 - 2 * r * cpsi + r**2
    m = 1 - r**2
    L = 1 - 3 * r * cpsi + 2 * r**2
    omega_psi = r * (cpsi - r) / D
    omega_r = spsi / D
    D_psi = 2 * r * spsi
    I_psi = 2 * qv**2 * (
        1
        + m * cpsi / (r * D)
        - m * spsi * D_psi / (r * D**2)
        + (3 * r**2 - 1) * omega_psi / r**2
    )
    I_psi_target = 4 * qv**2 * m * (cpsi - r) ** 2 / D**2

    m_r = -2 * r
    D_r = 2 * (r - cpsi)
    primitive_r = spsi * (
        (m_r * r * D - m * (D + r * D_r)) / (r * D) ** 2
    ) + 6 * omega / r + (3 * r**2 - 1) * (
        omega_r / r**2 - 2 * omega / r**3
    )
    I_r = 2 * qv**2 * primitive_r
    I_r_target = 4 * qv**2 * (
        omega / r**3 - m * spsi * L / (r**2 * D**2)
    )

    t = 2 * qv * sp.sqrt(m) * (cpsi - r) / D
    t_r = sp.diff(t, r)
    t_u_target = 2 * qv * m / D**2 * (
        m * ((cpsi - r) ** 2 - spsi**2)
        - r * (cpsi - r) * D
    )
    t_psi = -spsi * sp.diff(t, cpsi)
    t_psi_target = -2 * qv * m ** sp.Rational(3, 2) * spsi / D**2

    unit_circle = cpsi**2 + spsi**2 - 1

    def unit_circle_residual(expression: Any) -> str:
        numerator, _ = sp.together(sp.cancel(expression)).as_numer_denom()
        remainder = sp.rem(
            sp.Poly(numerator, spsi),
            sp.Poly(unit_circle, spsi),
        ).as_expr()
        return str(sp.factor(sp.cancel(remainder)))

    return {
        "domain": (
            "r != 0; m=1-r^2>0; D>0; 1-r*c>0; "
            "c^2+s^2=1; omega_psi=r*(c-r)/D; omega_r=s/D"
        ),
        "chain_rules": "c_psi=-s; s_psi=c; D_psi=2*r*s",
        "I_psi_residual_mod_unit_circle": unit_circle_residual(
            I_psi - I_psi_target
        ),
        "I_r_residual": str(sp.cancel(I_r - I_r_target)),
        "t_u_residual_mod_unit_circle": unit_circle_residual(
            t_r - t_u_target / m ** sp.Rational(3, 2)
        ),
        "t_psi_residual": str(sp.cancel(t_psi - t_psi_target)),
        "r_u": "m^(3/2)",
        "r0_limit": {
            "t": "2*q*cos(psi)",
            "t_u": "2*q*cos(2*psi)",
            "t_psi": "-2*q*sin(psi)",
            "I": "2*q^2*psi+q^2*sin(2*psi)",
            "I_u": "4*q^2*(sin(psi)+sin(3*psi)/3)",
        },
    }


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run fresh closed-form checks and return a serializable fixture.

    The module LF/EOF/AST preflight is deliberately first.  Numerical checks
    use only the closed formulas and fresh ``mpmath`` values; no saved source
    packet or producer is executed.  This function is opt-in and is never
    called during import.
    """

    preflight = _module_source_preflight(repo_root)
    source_binding = c5_loop_source_binding(repo_root)
    proof = _fresh_symbolic_proof()
    checks: dict[str, bool] = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "symbolic_I_psi_residual_zero": (
            proof["I_psi_residual_mod_unit_circle"] == "0"
        ),
        "symbolic_I_r_residual_zero": proof["I_r_residual"] == "0",
        "symbolic_t_u_residual_zero": (
            proof["t_u_residual_mod_unit_circle"] == "0"
        ),
        "symbolic_t_psi_residual_zero": (
            proof["t_psi_residual"] == "0"
        ),
        "symbolic_domain_and_r0_limit_present": (
            bool(proof.get("domain"))
            and proof.get("chain_rules") == "c_psi=-s; s_psi=c; D_psi=2*r*s"
            and isinstance(proof.get("r0_limit"), dict)
        ),
    }
    cases: dict[str, Any] = {}
    old_prec = mp.mp.prec

    def contains(interval: Any, value: Any) -> bool:
        return _lower(interval) <= value <= _upper(interval)

    try:
        c = MPIntervalContext()
        c.dps = 120

        # Keep this high-precision MPF regression under global 53-bit precision
        # while preserving the exact MPF tuple through the API boundary.
        with mp.workprec(240):
            high = mp.mpf(
                "0.123456789012345678901234567890123456789012345678901234567890123"
            )
        exact = mp.make_mpf(tuple(high._mpf_))
        mp.mp.prec = 53
        boxed = _coerce(c, high)
        with mp.workprec(240):
            checks["global_mp_prec53_high_precision_mpf_exact_tuple"] = (
                _lower(boxed) <= exact <= _upper(boxed)
                and abs(exact - _lower(boxed)) < mp.mpf("1e-100")
                and abs(_upper(boxed) - exact) < mp.mpf("1e-100")
            )
        mp.mp.prec = old_prec

        ordinary = closed_phase(c, Fraction(3, 4), 1, 2 * c.pi)
        checks["ordinary_upper_interval_not_collapsed"] = not (
            _lower(ordinary["Phi"]) == 1
            and _upper(ordinary["Phi"]) == 1
        )
        checks["ordinary_upper_phi_contains_one"] = contains(
            ordinary["Phi"], 1
        )
        checks["ordinary_upper_Iu_contains_zero"] = contains(
            ordinary["I_u"], 0
        )

        named = closed_phase(c, Fraction(3, 4), 1, "2*pi")
        checks["explicit_named_upper_normalized"] = (
            named["Phi"] == c.mpf(1)
            and named["I_u"] == c.mpf(0)
            and named["t_psi"] == c.mpf(0)
        )
        zero = closed_phase(c, Fraction(3, 4), 1, 0)
        checks["zero_singleton_normalized"] = (
            zero["Phi"] == c.mpf(0)
            and zero["I_u"] == c.mpf(0)
        )

        # Independent fresh closed-form reference.  This intentionally uses
        # no quadrature or saved source; ``mp.diff`` checks the primitive and
        # parameter derivatives against high-precision direct formulas.
        mp.mp.dps = 180
        qv = mp.mpf(1)
        uv = mp.mpf(3) / 4
        psiv = mp.mpf(11) / 10

        def reference_t(uu: mp.mpf, pp: mp.mpf) -> mp.mpf:
            ss = mp.sqrt(1 + uu * uu)
            rr = uu / ss
            return 2 * qv / ss * (mp.cos(pp) - rr) / (
                1 - 2 * rr * mp.cos(pp) + rr * rr
            )

        def reference_I(uu: mp.mpf, pp: mp.mpf) -> mp.mpf:
            ss = mp.sqrt(1 + uu * uu)
            rr = uu / ss
            mm = 1 - rr * rr
            dd = 1 - 2 * rr * mp.cos(pp) + rr * rr
            ww = mp.atan2(rr * mp.sin(pp), 1 - rr * mp.cos(pp))
            return 2 * qv * qv * (
                pp + mm * mp.sin(pp) / (rr * dd)
                + (3 * rr * rr - 1) * ww / (rr * rr)
            )

        reference = {
            "t": reference_t(uv, psiv),
            "t_u": mp.diff(lambda uu: reference_t(uu, psiv), uv),
            "t_psi": mp.diff(lambda pp: reference_t(uv, pp), psiv),
            "I": reference_I(uv, psiv),
            "I_u": mp.diff(lambda uu: reference_I(uu, psiv), uv),
        }
        out = closed_phase(c, Fraction(3, 4), 1, Fraction(11, 10))
        checks["primitive_and_t_derivatives_contained"] = all(
            contains(out[key], value) for key, value in reference.items()
        )
        K = 1 / (2 * mp.pi * (1 + 2 * qv * qv))
        contracted = (
            (1 + reference["t"] ** 2) * reference["t_u"]
            - reference["t_psi"] * reference["I_u"]
        )
        fixed = contracted / (1 + reference["t"] ** 2)
        normalized_b = (
            -mp.mpf(4) / 5 * (1 + 2 * qv * qv) / 2
            * (1 - reference["t"] ** 2)
            / (1 + reference["t"] ** 2) ** 3
            * contracted
        )
        phi_u = K * reference["I_u"]
        signed_left = (
            K * (1 + reference["t"] ** 2) * reference["t_u"]
            - reference["t_psi"] * phi_u
        )
        checks["lambda_tu_minus_tpsi_Phiu_equals_K_times_C"] = (
            abs(signed_left - K * contracted) < mp.mpf("1e-160")
        )
        checks["fixedphase_contraction_contained"] = contains(
            out["fixedphase_t_u"], fixed
        ) and contains(out["contractedC"], contracted)
        checks["bphi_normalized_derivative_contained"] = contains(
            out["bphi_u_over_e"], normalized_b
        )

        try:
            closed_phase(
                c,
                (Fraction(-1, 10), Fraction(1, 10)),
                1,
                Fraction(11, 10),
            )
        except ValueError as exc:
            checks["non_singleton_r_crossing_zero_rejected"] = (
                "contains zero" in str(exc)
            )
        else:
            checks["non_singleton_r_crossing_zero_rejected"] = False

        cases = {
            "u0_q1_psi_11_over_10": {
                "Phi": _interval_record(
                    closed_phase(c, 0, 1, Fraction(11, 10))["Phi"]
                ),
            },
            "u_3_over_4_q1_psi_11_over_10": {
                "Phi": _interval_record(out["Phi"]),
                "contractedC": _interval_record(out["contractedC"]),
            },
            "u_3_over_4_q1_psi_2pi_ordinary_iv": {
                "Phi": _interval_record(ordinary["Phi"]),
                "Phi_contains_one": checks[
                    "ordinary_upper_phi_contains_one"
                ],
                "Phi_is_not_singleton": checks[
                    "ordinary_upper_interval_not_collapsed"
                ],
            },
            "u_3_over_4_q1_psi_2pi_explicit_named": {
                "Phi": _interval_record(named["Phi"]),
                "I_u": _interval_record(named["I_u"]),
            },
        }
    finally:
        mp.mp.prec = old_prec

    passed = all(checks.values())
    return {
        "passed": passed,
        "module_sha256": preflight["sha256"],
        "source_binding": source_binding,
        "checks": checks,
        "cases": cases,
        "proof": proof,
        "scope": {
            "source_executed": False,
            "signed_source_actual_parameters_bound": (
                SIGNED_SOURCE_ACTUAL_PARAMETERS_BOUND
            ),
            "logcap_used_as_actual_value": False,
        },
    }


def fixture_bound_check(
    fixture: dict[str, Any],
    repo_root: str | Path | None = None,
) -> dict[str, Any]:
    """Check a previously saved fresh fixture against this final module."""

    preflight = _module_source_preflight(repo_root)
    source_binding = c5_loop_source_binding(repo_root)
    fixture_checks = fixture.get("checks") if isinstance(fixture, dict) else None
    checks = dict(fixture_checks) if isinstance(fixture_checks, dict) else {}
    checks["module_sha256_bound"] = (
        isinstance(fixture, dict)
        and fixture.get("module_sha256") == preflight["sha256"]
    )
    checks["source_binding_bound"] = (
        isinstance(fixture, dict)
        and fixture.get("source_binding") == source_binding
    )
    checks["fixture_passed"] = (
        isinstance(fixture, dict) and fixture.get("passed") is True
    )
    proof = fixture.get("proof") if isinstance(fixture, dict) else None
    checks["symbolic_proof_present"] = isinstance(proof, dict)
    if isinstance(proof, dict):
        checks["symbolic_proof_residuals_zero"] = all(
            proof.get(key) == "0"
            for key in (
                "I_psi_residual_mod_unit_circle",
                "I_r_residual",
                "t_u_residual_mod_unit_circle",
                "t_psi_residual",
            )
        )
    else:
        checks["symbolic_proof_residuals_zero"] = False
    return {
        "passed": all(checks.values()),
        "module_sha256": preflight["sha256"],
        "source_binding": source_binding,
        "checks": checks,
    }


__all__ = [
    "C5_LOOP_FUNCTION",
    "C5_LOOP_FUNCTION_AST_SHA256",
    "C5_LOOP_OUTPUT_PORTS",
    "C5_SOURCE",
    "C5_SOURCE_SHA256",
    "FIXTURE_RELATIVE_PATH",
    "MODULE",
    "SIGNED_SOURCE_ACTUAL_PARAMETERS_BOUND",
    "c5_loop_source_binding",
    "closed_phase",
    "fixture_bound_check",
    "fresh_fixture",
]
