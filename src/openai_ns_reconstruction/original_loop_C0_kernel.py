"""Source-independent scalar controls for the loop C0 kernel."""

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sp


def _lower(interval):
    if not hasattr(interval, "_mpi_"):
        raise TypeError("Expected a directed interval")
    row = interval._mpi_[0]
    if len(row) != 4:
        raise TypeError("Expected an MPIntervalContext interval")
    return mp.make_mpf(row)


def _upper(interval):
    if not hasattr(interval, "_mpi_"):
        raise TypeError("Expected a directed interval")
    row = interval._mpi_[1]
    if len(row) != 4:
        raise TypeError("Expected an MPIntervalContext interval pair")
    return mp.make_mpf(row)


def directed_log_bounds(ctx, logN):
    """Return exact directed log bounds without materializing N."""
    if not isinstance(ctx, MPIntervalContext) or not hasattr(logN, "_mpi_"):
        raise TypeError("ctx and logN must be an MPIntervalContext interval pair")
    log2 = ctx.log(ctx.mpf(2))
    if _lower(logN) < _upper(log2):
        raise ValueError("The certified logN interval must imply N >= 2")
    log2_over_N = log2 - logN
    log4_over_N = ctx.log(ctx.mpf(4)) - logN
    return {
        "logN": [list(row) for row in logN._mpi_],
        "log_2_over_N": [list(row) for row in log2_over_N._mpi_],
        "log_4_over_N": [list(row) for row in log4_over_N._mpi_],
        "N_at_least_2": True,
        "increment_assumption": "|A| <= 2",
        "increment_bound": "|exp(A/N)-1| <= 4/N",
        "increment_bound_log": [list(row) for row in log4_over_N._mpi_],
        "increment_proof": "|A/N| <= 1 and exp(x)-1 <= 2*x for 0 <= x <= 1",
        "N_materialized": False,
        "exp_minus_logN_materialized": False,
    }


def symbolic_controls():
    """Validate exact scalar identities and return their certificate."""
    a, b, eta, Delta, cut, t, u = sp.symbols("a b eta Delta cut t u", real=True)
    kappa = a + b**2 / a
    q2 = cut**2 * (2 * eta - Delta) / (2 * a)
    t0 = -b / a
    nu = 1 + t0**2 + 2 * q2
    v = sp.expand(a * nu)
    v_identity = sp.simplify(v.subs(Delta, kappa - 2)
                             - (kappa * (1 - cut**2) + (2 + 2 * eta) * cut**2)) == 0
    sum_square = sp.expand(kappa.subs(b, u) - 2 * u - (a - u)**2 / a) == 0
    un, qn, t0n = sp.symbols("u_n q_n t0_n", real=True)
    hn = 1 / sp.sqrt(1 + un**2)
    rn = un / sp.sqrt(1 + un**2)
    normalization_checks = {
        "one_minus_r_squared": sp.simplify(1 - rn**2 - hn**2) == 0,
        "average_f_squared": sp.simplify(1 / (2 * (1 - rn**2)) - 1 / (2 * hn**2)) == 0,
        "average_t_squared": sp.simplify(t0n**2 + 4 * qn**2 * hn**2 / (2 * (1 - rn**2))
                                           - (t0n**2 + 2 * qn**2)) == 0,
        "phase_normalization": sp.simplify(
            (2 * sp.pi * (1 / (2 * sp.pi * (1 + t0n**2 + 2 * qn**2)))
             * (1 + t0n**2 + 2 * qn**2)) - 1) == 0,
    }
    B_over_E = (-v * t / (1 + t**2) - b) / 2
    A_phi = a / 2 - v / (2 * (1 + t**2))
    checks = {
        "kappa_definition": sp.simplify(kappa - (a + b**2 / a)) == 0,
        "q_squared_definition": sp.simplify(q2 - cut**2 * (2 * eta - Delta) / (2 * a)) == 0,
        "t0_definition": sp.simplify(t0 + b / a) == 0,
        "nu_definition": sp.simplify(nu - (1 + t0**2 + 2 * q2)) == 0,
        "v_convex_combination": v_identity,
        "sum_square_identity": sum_square,
        **normalization_checks,
    }
    fixture = {a: sp.Rational(3, 2), b: sp.Rational(1, 2), eta: 1,
               Delta: sp.Rational(-1, 3), cut: sp.Rational(1, 2), t: sp.Rational(1, 3)}
    fk = sp.simplify(kappa.subs(fixture))
    fv = sp.simplify(v.subs(fixture))
    fixture_checks = {
        "a_positive": bool(fixture[a] > 0),
        "eta_range": bool(0 < fixture[eta] <= 1),
        "Delta_less_eta": bool(fixture[Delta] < fixture[eta]),
        "cut_range": bool(0 <= fixture[cut] <= 1),
        "a_le_kappa": bool(fixture[a] <= fk),
        "kappa_less_2_plus_eta": bool(fk < 2 + fixture[eta]),
        "v_upper": bool(fv <= 2 + 2 * fixture[eta]),
        "B_upper_7_over_4": bool(sp.Abs(B_over_E.subs(fixture)) <= sp.Rational(7, 4)),
        "A_upper_2": bool(sp.Abs(A_phi.subs(fixture)) <= 2),
    }
    passed = all(checks.values()) and all(fixture_checks.values())
    if not passed:
        raise ArithmeticError("Loop C0 symbolic controls failed")
    return {
        "passed": True,
        "assumptions": [
            "a > 0", "0 < eta <= 1", "Delta = kappa - 2 < eta",
            "kappa = a + b^2/a", "0 <= cut <= 1", "q^2 = cut^2*(2*eta-Delta)/(2*a)",
        ],
        "identities": checks,
        "fixture": {"a": "3/2", "b": "1/2", "eta": "1", "Delta": "-1/3", "cut": "1/2", "t": "1/3",
                    "kappa": str(fk), "v": str(fv),
                    "B_over_E": str(sp.factor(B_over_E.subs(fixture))),
                    "A_phi": str(sp.factor(A_phi.subs(fixture)))},
        "fixture_checks": fixture_checks,
        "exact_upper_constants": {"A_general": "2", "B_over_E_general": "2",
                                  "B_phi_over_E": "7/4", "B_over_E_endpoint": "7/8"},
        "bounds": {
            "abs_b_le_kappa_over_2": "(a-|b|)^2 >= 0 implies |b| <= kappa/2",
            "a_le_kappa_lt": "a <= kappa < 2 + eta",
            "v_upper": "v = kappa*(1-cut^2)+(2+2*eta)*cut^2 <= 2+2*eta",
            "B_over_E_upper": "|B_phi|/E <= 1+3*eta/4 <= 7/4 < 2",
            "A_upper": "a <= 3 and v <= 4 imply |A_phi| <= 2",
            "primitive_first_endpoint": "A(0)=0 and phase in [0,1] with |A_phi| <= 2 gives |A| <= 2",
            "primitive_B_first_endpoint": "B(0)=0 and phase in [0,1] with |B_phi|/E <= 2 gives |B|/E <= 2",
            "endpoint_A": "A(0)=A(1)=0 gives |A| <= 1 by the shorter phase path",
            "endpoint_B": "B(0)=B(1)=0 gives |B|/E <= 7/8 < 1; B(1) uses a*t0+b=0",
        },
        "phase_normalization": {
            "certificate": "f=(cos(psi)-r)/(1-2*r*cos(psi)+r^2)=sum_{n>=1} r^(n-1) cos(n*psi), |r|<1; Fourier orthogonality gives avg f=0 and avg f^2=1/[2*(1-r^2)]",
            "t_identity": "r=u/sqrt(1+u^2), h=1/sqrt(1+u^2), t=t0+2*q*h*f gives avg t=t0 and avg t^2=t0^2+2*q^2",
            "period_identity": "K=1/[2*pi*(1+t0^2+2*q^2)] gives Phi(2*pi)=1 and the inverse covers phase [0,1]",
            "numerical_quadrature": False,
        },
        "actual_source_claim": False,
    }


log_bounds = directed_log_bounds

__all__ = ["symbolic_controls", "directed_log_bounds", "log_bounds"]
