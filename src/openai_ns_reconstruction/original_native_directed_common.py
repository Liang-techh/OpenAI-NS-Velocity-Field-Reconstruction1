"""Changed-input directed formulas shared by the native-control repair.

Saved signed intervals are reused. Legacy logarithmic caps are not input
premises, except the explicitly repaired local inverse-L block. No accepted
source constructor, builder or value-jet recurrence is executed.
"""
from __future__ import annotations

import math
from mpmath import make_mpf
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, require
from .original_directed_abs_log import abs_log_upper
from .original_log_mixed_jets import _logsum
from .original_log_Z_phase_bounds import LogZ
from .original_native_log_jet4 import exp_relative, log_product, log_scale


def lower(value):
    return make_mpf(value._mpi_[0])


def upper(value):
    return make_mpf(value._mpi_[1])


def point(context, endpoint):
    return context.mpf([endpoint, endpoint])


def logs(context, values):
    return [abs_log_upper(context, value) for value in values]


def pack(rows):
    return [dict(ordinary_Z_order=k, exact_zero=value is None,
                 log_absolute_upper_bound=value) for k, value in enumerate(rows)]


def logsum(context, *values):
    return _logsum(context, [value for value in values if value is not None])


def maxlog(context, values):
    values = [value for value in values if value is not None]
    return None if not values else point(context, max(upper(value) for value in values))


def history(context, incoming, density, rate):
    """Positive Duhamel convex bound; density weight includes 1/rate."""
    require(type(rate) is int and rate > 0, "Positive source rate required")
    require(len(incoming) == len(density) == 5, "Five history rows required")
    log_rate = context.ln(context.mpf(rate))
    return [maxlog(context, (old, None if fresh is None else fresh-log_rate))
            for old, fresh in zip(incoming, density)]


def signed_rows(context, rows):
    require(isinstance(rows, list) and len(rows) == 5, "Five signed source rows required")
    require([row.get("ordinary_Z_order", k) for k, row in enumerate(rows)] == list(range(5)),
            "Signed source derivative orders changed")
    return [read(context, row.get("value", row)) for row in rows]


def cap_rows(context, rows):
    require(isinstance(rows, list) and len(rows) == 5, "Five repaired cap rows required")
    require([row["ordinary_Z_order"] for row in rows] == list(range(5)), "Cap orders changed")
    return [None if row["exact_zero"] else read(context, row["log_absolute_upper_bound"])
            for row in rows]


def value_product(context, left, right):
    """Explicit finite Leibniz formula, separate from the old value recurrence."""
    require(len(left) == len(right) == 5, "Five signed product rows required")
    return [sum((context.mpf(math.comb(k, j))*left[j]*right[k-j]
                 for j in range(k+1)), context.mpf(0)) for k in range(5)]


def analytic_log_jet(context, signed5):
    """Closed logarithmic derivatives of a positive signed five-row jet."""
    require(len(signed5) == 5 and lower(signed5[0]) > 0, "Positive signed value required")
    f0, f1, f2, f3, f4 = signed5
    r1, r2, r3, r4 = f1/f0, f2/f0, f3/f0, f4/f0
    return [context.ln(f0), r1, r2-r1*r1,
            r3-3*r1*r2+2*r1**3,
            r4-4*r1*r3-3*r2*r2+12*r1*r1*r2-6*r1**4]


def build_inputs(context, parent):
    """Bind corrected common factors to saved signed data, not old log caps."""
    require(isinstance(context, MPIntervalContext), "Interval context required")
    z = context.mpf(["37/100", "93/250"])
    delta = read(context, parent["actual_core_signed_direction"]["parameters"]["delta"])
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    log_f0 = read(context, amp["log_F0"])
    log_pstar = read(context, amp["log_Pstar"])
    log_lambda = read(context, parent["actual_macro_background_I1_numeric"]
                      ["geometric_measure"]["logLambda"])
    gradient = parent["actual_long_reference_native_controls"]["anchored_gradient"]
    require(gradient["ordinary_not_taylor"] is True
            and len(gradient["G_Z1_through_Z4"]) == 4, "Saved signed G1..G4 required")
    gradient_signed = [read(context, row) for row in gradient["G_Z1_through_Z4"]]
    gradient_logs = logs(context, gradient_signed)
    lambda_gradient = [None if value is None else value+log_lambda for value in gradient_logs]
    f2_ell = [None if value is None else value+context.ln(2) for value in lambda_gradient]
    f2_logs = log_scale(context, exp_relative(context, f2_ell),
                        2*point(context, upper(log_f0)))
    q = 1+z*z
    # These are ordinary derivatives of log(1+Z²), not Taylor coefficients.
    q_log_values = [context.ln(q), 2*z/q, 2*(1-z*z)/q**2,
                    4*z*(z*z-3)/q**3, -12*(z**4-6*z*z+1)/q**4]
    p0_source = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_pressure"]
    require(p0_source["ordinary_Z_orders"][:5] == list(range(5)), "Independent P0 Z0..4 required")
    p0_signed = [read(context, row) for row in p0_source["p0_jets"][:5]]
    p0_logs = [None if value is None else value+2*log_pstar for value in logs(context, p0_signed)]
    parameters = parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]["actual_logs"]
    local = parent["actual_native_control_directed_repair"]["charts"]
    inverse_l = cap_rows(context, local["restoration:post_restore"]
                        ["directed_L_block"]["directed_inverse_L_rows"])
    return dict(z=z, delta=delta, logLambda=log_lambda, logF0=log_f0,
                logPstar=log_pstar, logN=read(context, parameters["N"]),
                logeta=read(context, parameters["eta"]), logdstar=read(context, parameters["dstar"]),
                gradient_logs4=gradient_logs, lambda_gradient_logs4=lambda_gradient,
                gradient_signed4=gradient_signed, F2_logs5=f2_logs, P0_logs5=p0_logs,
                q_log_values5=q_log_values, q_derivative_logs4=logs(context, q_log_values[1:]),
                inverse_L_logs5=inverse_l,
                input_receipts=dict(saved_signed_gradient_reused=True,
                    old_gradient_recurrence_executed=False, unsafe_log_caps_used=False,
                    repaired_inverse_L_packet_used=True, source_handle_numbers_selected=False))


def p2(context, native, F2, P0, z, delta, logR, inverse_Utheta5, inverseL5):
    """Full native pressure formula with corrected independent factors."""
    alg = LogZ(context)
    jet = lambda rows: LogZ(context, list(rows[:4]))
    one = alg.constant(context.mpf(0))
    zz = jet(logs(context, [z, context.mpf(1), context.mpf(0), context.mpf(0)]))
    d = alg.add(one, alg.mul(zz, zz))
    h, ld, lopd = logs(context, [1-delta, delta, 1+delta])
    U, M, A = (jet(native[name]) for name in ("U", "M", "A"))
    Mz, Az = jet(native["M"][1:]), jet(native["A"][1:])
    hydro = alg.add(U,
        alg.scale(alg.add(M, alg.mul(zz, Mz)), h-context.ln(2)),
        alg.mul(U, alg.add(alg.scale(alg.mul(zz, M), h), alg.mul(d, Mz))),
        alg.scale(alg.mul(zz, A), context.ln(2)+ld), alg.mul(d, Az))
    pressure = alg.add(alg.scale(alg.mul(zz, jet(P0)), context.ln(2)+lopd),
                       alg.mul(d, jet(P0[1:])))
    fb, fc = log_product(context, F2, native["B"]), log_product(context, F2, native["C"])
    dressed = alg.add(alg.scale(alg.mul(zz, jet(fb)), context.ln(2)+ld),
        alg.mul(d, jet(fb[1:])), alg.scale(alg.mul(zz, jet(fc)), context.ln(2)+lopd),
        alg.mul(d, jet(fc[1:])))
    numerator = alg.add(hydro, pressure, alg.scale(dressed, logR))
    result = alg.scale(alg.mul(alg.mul(numerator, jet(inverse_Utheta5)), jet(inverseL5)), logR)
    values = lambda actual: [actual[k] for k in range(4)]
    return values(result), dict(hydro=values(hydro), independent_pressure=values(pressure),
        dressed_swirl=values(dressed), numerator=values(numerator),
        inverse_L=list(inverseL5[:4]), inverse_Utheta=list(inverse_Utheta5[:4]),
        formula="p2=R*Nnative/(Utheta*L)", full_P0_and_P0_Z=True,
        all_M_ABC_terms_retained=True, physical_J_not_in_denominator=True,
        signed_value_selected=False)


def fresh_fixture():
    """New exact rational derivatives and independent pressure sectors."""
    import sympy as s
    context = MPIntervalContext()
    context.prec = 1600
    fixtures = (([1, 1, 0, 0, 0], [0, 1, -1, 2, -6]),
                ([1, 0, 2, 0, 0], [0, 0, 2, 0, -12]))
    for values, expected in fixtures:
        actual = analytic_log_jet(context, [context.mpf(v) for v in values])
        require(all(lower(v) <= e <= upper(v) for v, e in zip(actual, expected)),
                "Closed log-derivative fixture failed")
    z = s.Symbol("Z")
    location, de, radius = s.Rational(1, 3), s.Rational(1, 7), s.Integer(3)
    q, L = 1+z*z, 1-de*z*z

    def scalar(value):
        value = s.Rational(value)
        return context.mpf(int(value.p))/int(value.q)

    def jet(function):
        return [scalar(s.diff(function, z, k).subs(z, location)) for k in range(5)]

    inv_u, inv_l = logs(context, jet(q/5)), logs(context, jet(1/L))
    sectors = (
        dict(U=0, M=0, A=0, B=0, C=0, F=0, P=1+z+z**4),
        dict(U=0, M=0, A=1+z**3, B=0, C=0, F=0, P=0),
        dict(U=0, M=0, A=0, B=1+z**2, C=0, F=(1+z)**2, P=0),
        dict(U=0, M=0, A=0, B=0, C=1+z**3, F=(1+z)**2, P=0),
        dict(U=0, M=1+z**2, A=0, B=0, C=0, F=0, P=0),
        dict(U=1+z, M=2+z*z, A=3+z**3, B=1+z*z, C=2+z**3,
             F=(1+z)**2, P=1+z**4),
    )
    nonzero = 0
    for sector in sectors:
        U, M, A, B, C, F, P = (s.sympify(sector[name]) for name in ("U", "M", "A", "B", "C", "F", "P"))
        d, h = 1-z*z, 1-de
        # Independently differentiate the signed rational function, rather
        # than the logarithmic majorant graph used by p2().
        numerator = (-U+h*(M-z*s.diff(M, z))/2
            +U*(h*z*M+d*s.diff(M, z))+2*de*z*A-d*s.diff(A, z)
            +2*(1+de)*z*P-d*s.diff(P, z)
            +radius*(-2*de*z*F*B+d*s.diff(F*B, z)
                      +2*(1+de)*z*F*C-d*s.diff(F*C, z)))
        actual_function = radius*numerator*q/(5*L)
        native = {name: logs(context, jet(s.sympify(sector[name])))
                  for name in ("U", "M", "A", "B", "C")}
        bound, _ = p2(context, native, logs(context, jet(F)), logs(context, jet(P)),
                      scalar(location), scalar(de), context.ln(3), inv_u, inv_l)
        for k, cap in enumerate(bound):
            expected = s.cancel(s.diff(actual_function, z, k).subs(z, location))
            if expected == 0:
                continue
            reference = abs_log_upper(context, scalar(expected))
            require(cap is not None and upper(cap) >= upper(reference),
                    "Independent signed pressure sector escaped its majorant")
            nonzero += 1
    return dict(passed=True, closed_log_fixtures=2, historical_inputs_used=False,
                independent_pressure_sectors=len(sectors), nonzero_pressure_rows=nonzero,
                accepted_self_checks_called=False)
