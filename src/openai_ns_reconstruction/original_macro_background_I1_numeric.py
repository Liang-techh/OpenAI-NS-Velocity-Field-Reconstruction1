"""Pre-finite-N native macro signed six-I1 C1 enclosures.

This module is a small arithmetic adapter over the saved current source
packets.  It consumes the C6 comparison fields and the saved whole-micro
native background, then propagates the latter over the macro chart with the
new ordinary jet bounds.  The finite-N graph source and its outer ``R_Q_y``
ports remain receipts for a later correction; they are deliberately not
evaluated here.

The macro coordinate is ``m`` in ``[0, 1]`` and

    R = Ra*exp(2*h_b + ell*m),  ell = Y - 2*h_b,
    Y = log(100) - log(Ra).

The saved width is ``loghb``.  We retain that logarithm and only use the
bound ``0 < h_b <= 2**-1024`` in interval arithmetic, so the astronomical
physical width, the finite-N scale, and F0 are never materialized.  The
dimensionless exponentials ``exp(-2*ell)``, ``exp(logLambda)``, and
``exp(2*logPstar)`` have representable integer-sized MPF exponents and are
evaluated where the reduced formulas require them.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

from .original_core_signed_direction import require, sha, read, record, source_guard
from .original_ordinary_jet_bounds import (
    jet_add,
    jet_div,
    jet_exp,
    jet_log,
    jet_mul,
    jet_neg,
    jet_reciprocal,
)


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_macro_correlated_m15.json"
PARENT_SHA256 = "8354af18754496450fb41ff0335e2c82f503c0d64852496da11ee1546fb23577"
B2_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_b_p2_mixed_bounds.json"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
CORE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
CORE_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json"
JET_SOURCE = "src/openai_ns_reconstruction/original_ordinary_jet_bounds.py"
SIGNED_HELPER = "src/openai_ns_reconstruction/original_core_signed_direction.py"

NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
HISTORY_NAMES = ("H", "K", "M", "A", "B", "C")
HISTORY_RATES = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
ORDER = 3


def _load(path):
    return json.loads(Path(path).read_bytes())


def _symmetric(c, value):
    radius = c.fabs(value).b
    require(radius >= 0, "A nonnegative interval radius is required")
    return c.mpf([-radius, radius])


def _rows(values):
    return [dict(ordinary_Z_order=k, value=record(value))
            for k, value in enumerate(values)]


def _row_values(c, rows, order=ORDER):
    require(isinstance(rows, list) and len(rows) >= order + 1,
            "Saved ordinary jet rows are incomplete")
    return [read(c, row.get("value", row)) for row in rows[:order + 1]]


def _pin(root, parent, path):
    """Pin every file read by this adapter through the supplied parent."""
    if path.startswith(BASE):
        expected = parent["input_hashes"].get(Path(path).name)
    else:
        expected = parent["implementation_hashes"].get(path)
    require(expected is not None, "Latest parent has no hash for " + path)
    require(sha(root / path) == expected, "Changed saved input: " + path)
    return dict(path=path, expected_sha256=expected, passed=True)


def _hull(c, left, right):
    """Interval hull of two ordinary-jet row enclosures."""
    return c.mpf([min(left.a, right.a), max(left.b, right.b)])


def _jet_scale(c, jet, scalar):
    return [scalar * value for value in jet]


def _jet_sub(c, left, right):
    return jet_add(c, left, jet_neg(c, right))


def _jet_const(c, value, order=ORDER):
    return [c.mpf(value)] + [c.mpf(0)] * order


def _jet_add_radius(c, jet, radius):
    return [value + _symmetric(c, radius) for value in jet]


def _saved_delta(c, parent):
    return read(c, parent["actual_core_signed_direction"]["parameters"]["delta"])


def _saved_logp(c, parent):
    return read(c, parent["actual_core_signed_direction"]["parameters"]["logPstar"])


def _source_guards(root):
    """AST-bind the complete C6 directions/fields and the core F0 recipe."""
    return dict(
        c6_densities=source_guard(root, C6_SOURCE, "densities", (
            "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",
        )),
        c6_histories=source_guard(root, C6_SOURCE, "histories", (
            "for name,rate in RATES.items():",
            "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,",
            "original_incoming_memory_retained=True,source_measure_applied_once=True",
        )),
        c6_comparison_directions=source_guard(root, C6_SOURCE, "comparison_directions", (
            "f,h,k=(b.mul(F0,jet) for jet in (phi,J['H'],J['K']))",
            "W=a.sub(one,a.add(a.mul(a.sub(one,de),a.mul(z,M)),a.mul(d,Mz)))",
            "D=a.div(a.add(a.neg(W),a.div(angular,a.scale(a.trim(f),2),",
        )),
        c6_actual_fields=source_guard(root, C6_SOURCE, "actual_fields", (
            "full_drive=a.add(a.scale(known['drive_hydro'],radius)",
            "Vy=a.scale(a.mul(ratio,full_drive),g.neg(chi))",
            "actual_prescribed_axial_equation=True,all_hydro_pressure_swirl_sectors_retained=True",
        )),
        core_F0_recipe=source_guard(root, CORE_SOURCE, "__init__", (
            "gradient=a.div(a.mul(L,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))",
            "G=controls.integral(g,gv,g.nodes[v.node]['name'],self.parameters['anchor'],z,",
            "Gjet=c3.C3Function(G,gradient.value,gradient.Z,gradient.ZZ)",
            "logLambda=g.add(g.mul(c(4),logP),c(1000))",
            "Lambda,epsilon=g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda))",
            "logF=a.add(a.fixed(g.neg(self.parameters['logC'])),a.neg(a.scale(Gjet,Lambda)))",
            "F0=diff.lift(source.C1Function(F0value,g.mul(F0value,logF.Z)))",
        )),
        c6_endpoint_inclusion=source_guard(root, C6_SOURCE, "__init__", (
            "self.comparison_micro_end=dict(phi=a.at(barphi,s,c(2)),U=a.at(barU,s,c(2)),",
            "barphi2,barU2=(self.comparison_micro_end[name] for name in ('phi','U'))",
            "actualmacro=actual_fields(b,macroknown,self.actual_micro_end,m,length,hb,Rmacro,F0,Pstar,barphi2,'macro')",
            "actualmacroJ=histories(b,actualmacro['phi'],actualmacro['U'],self.actual_moment_end,m,length,'actual_macro')",
        )),
    )


def _source_functional_identities():
    """Fresh symbolic checks for the C6 macro reductions used below."""
    import sympy as sp

    z, de, lam, phi = sp.symbols("z delta Lambda phi", positive=True)
    H, Hz, K, Kz, M, Mz, g = sp.symbols("H H_Z K K_Z M M_Z g")
    U, A, Az, B, Bz, C, Cz = sp.symbols("U A A_Z B B_Z C C_Z")
    P, Pz, R, F2, F2z = sp.symbols("P_0 P_0_Z R F_0^2 (F_0^2)_Z", positive=True)
    d = 1 - z**2
    L = 1 - de*z**2
    h = 1 - de
    Ared = (1 - de / 2) * H - h*z*Hz/2 - Kz + (2*de - 1)*z*K
    Bred = h*z*H/2 + d*K
    Wred = 1 - h*z*M - d*Mz
    Dred = (-Wred + (Ared + lam*g*Bred)/(2*phi)) / L
    hydro = (-Wred*U + h*(M - z*Mz)/2 + 2*de*z*A - d*Az) / (2*L)
    pressure = (2*(1 + de)*z*P - d*Pz) / (2*L)
    F2Bz = F2z*B + F2*Bz
    F2Cz = F2z*C + F2*Cz
    dressed = (-2*de*z*(F2*B) + d*F2Bz
                + 2*(1 + de)*z*(F2*C) - d*F2Cz) / (2*L)
    quotient = sp.cancel(R**2 * F2 * (dressed / F2) - R**2 * dressed)
    q = (2*z*U - h*z*M - d*Mz) / L
    checks = dict(
        A_reduction=sp.simplify(Ared - ((1 - de / 2)*H - h*z*Hz/2 - Kz + (2*de - 1)*z*K)) == 0,
        B_reduction=sp.simplify(Bred - (h*z*H/2 + d*K)) == 0,
        W_reduction=sp.simplify(Wred - (1 - h*z*M - d*Mz)) == 0,
        D_over_R_reduction=sp.simplify(Dred*L - (-Wred + (Ared + lam*g*Bred)/(2*phi))) == 0,
        Q_recovery=sp.simplify(q*L - (2*z*U - h*z*M - d*Mz)) == 0,
        F0_squared_external_internal_cancellation=quotient == 0,
    )
    require(all(checks.values()), "New macro source identity failed")
    return dict(
        passed=True,
        checks=checks,
        ordinary_derivatives_have_no_factorials=True,
        C6_comparison_reduction=(
            "A=(1-delta/2)H-(1-delta)Z*H_Z/2-d*K_Z+(2*delta-1)Z*K; "
            "B=(1-delta)Z*H/2+d*K; W=1-(1-delta)Z*M-d*M_Z"
        ),
        D_over_R="(-W+(A+Lambda*gaxis*B)/(2*phi))/L",
        full_drive=(
            "R*hydro+R*Pstar^2*pressure+R^2*dressed_swirl"
        ),
        dressed_swirl_product_derivatives=(
            "(F0^2*B)_Z=F0^2_Z*B+F0^2*B_Z; "
            "(F0^2*C)_Z=F0^2_Z*C+F0^2*C_Z"
        ),
        F0_log_derivative="F0_Z/F0=-Lambda*gaxis",
        F0_squared_cancellation_exact=True,
        Q_recovery="Q=(2*Z*U-(1-delta)*Z*M-d*M_Z)/L",
    )


def _comparison_histories(c, parent):
    """Build macro comparison history hulls from saved incoming rows."""
    ea = parent["actual_micro_Ea_mixed_bounds"]["comparison_source_Z0_Z4_enclosures"]
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    source_order = ORDER + 1
    phi = _row_values(c, ea["comparison_Phi"], source_order)
    U = _row_values(c, ea["comparison_U"], source_order)
    incoming = {
        "H": _row_values(c, ea["JH"], source_order),
        "K": _row_values(c, ea["JK"], source_order),
        "M": _row_values(c, ea["M"], source_order),
    }
    abc = b2["comparison_complete_ABC_history_rows"]
    for name in ("A", "B", "C"):
        incoming[name] = _row_values(c, abc[name], source_order)

    # The requested macro equilibrium normalization is the C6 comparison
    # normalization: H=phi, K=phi*U, M=U, A=U^2, B=phi^2/2, C=phi^2.
    phi_u = jet_mul(c, phi, U)
    phi2 = jet_mul(c, phi, phi)
    u2 = jet_mul(c, U, U)
    equilibrium = {
        "H": phi,
        "K": phi_u,
        "M": U,
        "A": u2,
        "B": [value / 2 for value in phi2],
        "C": phi2,
    }
    hulls = {}
    receipts = {}
    for name in HISTORY_NAMES:
        hull = [_hull(c, incoming[name][k], equilibrium[name][k])
                for k in range(source_order + 1)]
        hulls[name] = hull
        receipts[name] = dict(
            ordinary_Z_orders=list(range(source_order + 1)),
            rate=HISTORY_RATES[name],
            incoming=_rows(incoming[name]),
            equilibrium=_rows(equilibrium[name]),
            hull=_rows(hull),
            formula=(
                "equilibrium+exp(-rate*ell*m)*(incoming-equilibrium)"
            ),
            decay_interval="[0,1] for m in [0,1]",
            decay_Z_independent=True,
            convex_hull_enclosure=True,
            incoming_memory_retained=True,
            comparison_only=True,
        )
    return dict(phi=phi, U=U, incoming=incoming, equilibrium=equilibrium,
                hulls=hulls, receipts=receipts,
                source="actual_micro_Ea_mixed_bounds.comparison_source_Z0_Z4_enclosures",
                abc_source=("actual_micro_b_p2_mixed_bounds."
                             "comparison_complete_ABC_history_rows"))


def _macro_geometry(c, parent):
    """Return the Z-independent macro geometry and positive measure factor."""
    logp = _saved_logp(c, parent)
    loglambda = 4 * logp + 1000
    logra = c.log(c.mpf(4)) - loglambda
    Y = c.log(c.mpf(100)) - logra
    loghb = read(c, parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]["hb"])
    hb_threshold = -1024 * c.log(c.mpf(2))
    require(loghb.b < hb_threshold.a,
            "Saved loghb must be below log(2^-1024) before using the width cap")
    hb_upper = c.mpf(2) ** (-1024)
    # 0 < h_b <= 2^-1024; ell = Y - 2 h_b.
    ell = c.mpf([Y.a - 2 * hb_upper.b, Y.b])
    require(ell.a > 0, "Macro length must be positive")
    K = (1 - c.exp(-2 * ell)) / 2
    require(K.a > 0 and K.b <= c.mpf(1), "Positive macro geometry factor required")
    return dict(
        Wmacro=10000,
        Wmacro_record=record(c.mpf(10000)),
        logPstar=record(logp),
        logLambda=record(loglambda),
        logRa=record(logra),
        Y=record(Y),
        loghb=record(loghb),
        loghb_threshold=record(hb_threshold),
        saved_loghb_below_threshold=True,
        saved_hb_field="loghb",
        physical_hb_definition="h_b=exp(loghb)",
        hb_upper=record(hb_upper),
        ell=record(ell),
        K=record(K),
        R_bound=record(c.mpf([0, 100])),
        R_positive_for_0_lt_m_le_1=True,
        R_upper=100,
        radius_definition="R=Ra*exp(2*h_b+ell*m)",
        length_definition="ell=Y-2*h_b; Y=log(100)-log(Ra)",
        measure_formula=(
            "integral_0^1 R^2*ell*dm/Wmacro=(1-exp(-2*ell))/2"
        ),
        physical_R_squared_length_applied_once=True,
        radius_Z_independent=True,
        length_Z_independent=True,
        no_physical_hb_materialized=True,
        no_physical_Ra_materialized=True,
    )


def _comparison_directions(c, comparison, parent, geometry):
    """Recover C6 macro D_over_R with F0 eliminated before boxing."""
    ea = parent["actual_micro_Ea_mixed_bounds"]["comparison_source_Z0_Z4_enclosures"]
    phi = comparison["phi"][:ORDER + 1]
    U = comparison["U"][:ORDER + 1]
    H_all, K_all, M_all = (comparison["hulls"][name]
                           for name in ("H", "K", "M"))
    A_all, B_all, C_all = (comparison["hulls"][name]
                           for name in ("A", "B", "C"))
    H, K, M = H_all[:ORDER + 1], K_all[:ORDER + 1], M_all[:ORDER + 1]
    A, B, C = A_all[:ORDER + 1], B_all[:ORDER + 1], C_all[:ORDER + 1]
    z_saved = _row_values(c, ea["Z"], ORDER + 1)
    z = z_saved[:ORDER + 1]
    delta = _saved_delta(c, parent)
    one = _jet_const(c, 1)
    de = _jet_const(c, delta)
    d = _jet_sub(c, one, jet_mul(c, z, z))
    L = _jet_sub(c, one, _jet_scale(c, jet_mul(c, z, z), delta))
    require(L[0].a > 0, "Macro L=1-delta*Z^2 must be positive")
    # Consume the saved C0 comparison d/L rows as containment receipts.
    for key, current in (("d", d), ("L", L)):
        saved = _row_values(c, ea[key])
        require(all(saved[k].a <= current[k].a and current[k].b <= saved[k].b
                    for k in range(ORDER + 1)),
                "Saved comparison " + key + " rows do not contain the new algebra")

    h = _jet_const(c, 1 - delta)
    half = c.mpf("0.5")
    H_z, K_z, M_z = H_all[1:ORDER + 2], K_all[1:ORDER + 2], M_all[1:ORDER + 2]
    Ared = _jet_sub(c, _jet_scale(c, H, 1 - delta / 2),
                    _jet_scale(c, jet_mul(c, h, jet_mul(c, z, H_z)), half))
    Ared = _jet_sub(c, Ared, jet_mul(c, d, K_z))
    Ared = jet_add(c, Ared, _jet_scale(c, jet_mul(c, z, K), 2 * delta - 1))
    Bred = jet_add(c, _jet_scale(c, jet_mul(c, h, jet_mul(c, z, H)), half),
                   jet_mul(c, d, K))
    Wred = _jet_sub(c, one, jet_add(c, jet_mul(c, h, jet_mul(c, z, M)),
                                  jet_mul(c, d, M_z)))
    loglambda = read(c, geometry["logLambda"])
    Lambda = c.exp(loglambda)
    # Keep the Lambda factor Z-independent and explicit; inverse F0 is never
    # formed.  The compact expression is exactly the C6 F0 quotient.
    gaxis = read_jet(c, ea["gaxis"], ORDER + 1)[:ORDER + 1]
    lam_g = _jet_scale(c, gaxis, Lambda)
    angular = jet_div(c, jet_add(c, Ared, jet_mul(c, lam_g, Bred)),
                      _jet_scale(c, phi, 2))
    numerator = jet_add(c, jet_neg(c, Wred), angular)
    D = jet_div(c, numerator, L)
    return dict(
        Z=z,
        delta=delta,
        h=1 - delta,
        d=d,
        L=L,
        Lambda=Lambda,
        phi=phi,
        U=U,
        H=H,
        K=K,
        M=M,
        A=A,
        B=B,
        C=C,
        A_reduced=Ared,
        B_reduced=Bred,
        W_reduced=Wred,
        gaxis=gaxis,
        D_over_R=D,
        saved_d_rows=_rows(_row_values(c, ea["d"], ORDER + 1)),
        saved_L_rows=_rows(_row_values(c, ea["L"], ORDER + 1)),
        F0_inverse_not_materialized=True,
        F0_Z_over_F0="-Lambda*gaxis",
        ordinary_Z_orders=list(range(ORDER + 1)),
    )


def read_jet(c, rows, order=ORDER):
    """Read a saved record jet without importing any accepted numeric builder."""
    return _row_values(c, rows, order)


def _jet_neg(c, jet):
    return jet_neg(c, jet)


def _full_signed_axial_drive(c, comparison, directions, parent, geometry):
    """C6 actual_fields drive, with F0^2 kept as the saved coarse cap."""
    delta = directions["delta"]
    z, d, L = directions["Z"], directions["d"], directions["L"]
    U = comparison["U"][:ORDER + 1]
    M_all = comparison["hulls"]["M"]
    A_all, B_all, C_all = (comparison["hulls"][name]
                           for name in ("A", "B", "C"))
    M = M_all[:ORDER + 1]
    A, B, C = A_all[:ORDER + 1], B_all[:ORDER + 1], C_all[:ORDER + 1]
    M_z, A_z, B_z, C_z = (M_all[1:ORDER + 2], A_all[1:ORDER + 2],
                          B_all[1:ORDER + 2], C_all[1:ORDER + 2])
    W = directions["W_reduced"]
    one = _jet_const(c, 1)
    half = c.mpf("0.5")
    den = _jet_scale(c, L, 2)
    h = 1 - delta
    hydro_num = jet_add(c, jet_neg(c, jet_mul(c, W, U)),
                        _jet_scale(c, jet_mul(c, _jet_sub(c, M, jet_mul(c, z, M_z)),
                                             _jet_const(c, h)), half))
    hydro_num = jet_add(c, hydro_num, _jet_scale(c, jet_mul(c, z, A), 2 * delta))
    hydro_num = _jet_sub(c, hydro_num, jet_mul(c, d, A_z))
    hydro = jet_div(c, hydro_num, den)
    p0_all = _row_values(c, parent["actual_micro_b_p2_mixed_bounds"]
                         ["whole_tube_pressure"]["p0_jets"], ORDER + 1)
    p0 = p0_all[:ORDER + 1]
    pressure_num = _jet_scale(c, jet_mul(c, z, p0), 2 * (1 + delta))
    pressure_num = _jet_sub(c, pressure_num,
                            jet_mul(c, d, p0_all[1:ORDER + 2]))
    pressure = jet_div(c, pressure_num, den)
    f2 = [c.mpf([0, 1])] + [c.mpf([-1, 1])] * (ORDER + 1)
    f2_b_all = jet_mul(c, f2, B_all)
    f2_c_all = jet_mul(c, f2, C_all)
    f2_b, f2_c = f2_b_all[:ORDER + 1], f2_c_all[:ORDER + 1]
    dressed_num = _jet_scale(c, jet_mul(c, z, f2_b), -2 * delta)
    dressed_num = jet_add(c, dressed_num,
                          jet_mul(c, d, f2_b_all[1:ORDER + 2]))
    dressed_num = jet_add(c, dressed_num,
                          _jet_scale(c, jet_mul(c, z, f2_c), 2 * (1 + delta)))
    dressed_num = _jet_sub(c, dressed_num,
                           jet_mul(c, d, f2_c_all[1:ORDER + 2]))
    dressed_swirl = jet_div(c, dressed_num, den)
    R = [c.mpf([0, 100])] + [c.mpf(0)] * ORDER
    R2 = jet_mul(c, R, R)
    P2 = c.exp(2 * _saved_logp(c, parent))
    full = jet_add(c, _jet_scale(c, hydro, R[0]),
                   _jet_scale(c, _jet_scale(c, pressure, P2), R[0]))
    full = jet_add(c, full, jet_mul(c, R2, dressed_swirl))
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    require(amp["exact_log_F0"] == "-logC-Lambda*G; 0<=G<=Gbar on this tube" and
            amp["F0_squared_ordinary_Z0_Z4_coarse_absolute_cap_one"] is True and
            amp["selected_logC_not_replaced_by_minimum_guard"] is True,
            "Saved exact log F0 and coarse F0^2 receipt are required")
    return dict(
        hydro=hydro,
        pressure=pressure,
        dressed_swirl=dressed_swirl,
        full_signed_axial_drive=full,
        p0=p0,
        Pstar_squared=P2,
        F0_squared_cap=f2,
        exact_log_F0=amp["exact_log_F0"],
        log_F0=record(read(c, amp["log_F0"])),
        selected_logC_not_replaced_by_minimum_guard=True,
        F0_squared_coarse_cap_consumed=True,
        R_bound=R,
        R_squared_bound=R2,
        R_positive_upper=100,
        R_Z_independent=True,
        F0_squared_external_internal_cancellation=True,
        F0_inverse_not_materialized=True,
        source_units="p0=P0/Pstar^2; Pstar^2 applied exactly once in pressure drive",
        ordinary_Z_orders=list(range(ORDER + 1)),
    )


def _log_bound(c, loghb, ell, cap, divisor=2):
    require(cap > 0, "Positive product cap required for logarithmic error proof")
    return loghb + c.log(ell) + c.log(c.mpf(100)) + c.log(cap) - c.log(c.mpf(divisor))


def _actual_macro_fields(c, comparison, directions, drive, parent, geometry):
    """Propagate pre-finite-N native phi/U/M enclosures over macro m."""
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    native = b2["native_angular_and_axial_provider"]
    native_phi_all = _row_values(c, native["native_phi"], ORDER + 1)
    native_phi = native_phi_all[:ORDER + 1]
    numeric = parent["actual_micro_background_I1_numeric"]
    native_U_rows = native["native_U"][:ORDER + 1]
    saved_U_rows = numeric["signed_source_rows"]["U"][:ORDER + 1]
    require(native_U_rows == [row.get("value", row) for row in saved_U_rows],
            "B2 native_U must equal the saved signed whole-micro U rows")
    initial_U = _row_values(c, native_U_rows)
    initial_M = _row_values(c, numeric["signed_source_rows"]["M"])
    require(native_phi[0].a > 0, "Native macro phi initial lower bound must be positive")
    loghb = read(c, geometry["loghb"])
    ell = read(c, geometry["ell"])
    threshold = -1024 * c.log(c.mpf(2))
    D = directions["D_over_R"]
    phi_log_caps = [c.fabs(value).b for value in D]
    phi_log_bounds = [_log_bound(c, loghb, ell, cap, divisor=2)
                      for cap in phi_log_caps]
    require(all(bound.b < threshold.a for bound in phi_log_bounds),
            "Macro angular log error does not fit the saved 2^-1024 budget")
    radius = c.mpf(2) ** (-1024)
    initial_logphi = jet_log(c, native_phi)
    actual_logphi = _jet_add_radius(c, initial_logphi, radius)
    phi = jet_exp(c, actual_logphi)
    ratio = jet_div(c, phi, comparison["phi"][:ORDER + 1])
    transported = jet_mul(c, ratio, drive["full_signed_axial_drive"])
    transport_caps = [c.fabs(value).b for value in transported]
    transport_log_bounds = [_log_bound(c, loghb, ell, cap, divisor=1)
                            for cap in transport_caps]
    require(all(bound.b < threshold.a for bound in transport_log_bounds),
            "Macro axial transport error does not fit the saved 2^-1024 budget")
    U = _jet_add_radius(c, initial_U, radius)
    M = [_hull(c, initial_M[k], U[k]) for k in range(ORDER + 1)]
    return dict(
        initial_phi=native_phi,
        initial_U=initial_U,
        initial_M=initial_M,
        initial_U_source=(
            "actual_micro_b_p2_mixed_bounds.native_angular_and_axial_provider.native_U"
            "; exact-equal to actual_micro_background_I1_numeric.signed_source_rows.U"
        ),
        initial_M_source=(
            "actual_micro_background_I1_numeric.signed_source_rows.M"
            "; saved signed whole-micro endpoint source"
        ),
        initial_logphi=initial_logphi,
        actual_logphi=actual_logphi,
        phi=phi,
        U=U,
        M=M,
        phi_over_comparison=ratio,
        transported_product=transported,
        error_radius=radius,
        native_M_nonreset=True,
        native_M_not_identified_with_comparison_M=True,
        phi_positive_floor=record(phi[0].a),
        ordinary_Z_orders=list(range(ORDER + 1)),
        error_log_proofs=dict(
            threshold=record(threshold),
            phi_log_rows=[dict(ordinary_Z_order=k, maxabs_D=record(c.mpf([0, cap])),
                                log_upper=record(bound), passed=bound.b < threshold.a)
                          for k, (cap, bound) in enumerate(zip(phi_log_caps, phi_log_bounds))],
            U_transport_rows=[dict(ordinary_Z_order=k,
                                   maxabs_product=record(c.mpf([0, cap])),
                                   log_upper=record(bound), passed=bound.b < threshold.a)
                              for k, (cap, bound) in enumerate(zip(transport_caps, transport_log_bounds))],
            common_radius=record(radius),
            phi_formula="log(phi)-log(phi_initial)=-h_b/2*integral(R*D_over_R*dy)",
            U_formula="U_y=-h_b*(phi/barphi)*full_signed_axial_drive",
            dy_formula="dy=ell*dm",
            h_b_definition="exp(loghb)",
            no_hb_materialized=True,
        ),
    )


def _recover_Q(c, fields, directions):
    """Recover native pre-finite-N Q_b from the same linear C6 identity."""
    z = directions["Z"][:3]
    d = directions["d"][:3]
    L = directions["L"][:3]
    delta = directions["delta"]
    U, M = fields["U"][:3], fields["M"][:3]
    Mz = fields["M"][1:4]
    numerator = _jet_scale(c, jet_mul(c, z, U), 2)
    numerator = _jet_sub(c, numerator,
                         _jet_scale(c, jet_mul(c, z, M), 1 - delta))
    numerator = _jet_sub(c, numerator, jet_mul(c, d, Mz))
    Q = jet_div(c, numerator, L)
    require(len(Q) == ORDER, "Q recovery requires ordinary Z0..Z2")
    return Q


def _six_integrands(c, fields, Q):
    U = fields["U"]
    return dict(
        q=Q[:2],
        q_Z=Q[1:3],
        q_squared=jet_mul(c, Q[:2], Q[:2]),
        U_q_Z=jet_mul(c, U[:2], Q[1:3]),
        U_squared=jet_mul(c, U[:2], U[:2]),
        U_U_Z=jet_mul(c, U[:2], U[1:3]),
    )


def _macro_I1_and_D(c, integrands, directions, geometry):
    K = read(c, geometry["K"])
    macro_I1 = {name: [K * value for value in jet]
                for name, jet in integrands.items()}
    z = directions["Z"][:2]
    delta = directions["delta"]
    d = directions["d"][:2]
    L = directions["L"][:2]
    h = 1 - delta
    invL = jet_reciprocal(c, L)
    invL2 = jet_mul(c, invL, invL)
    half = [c.mpf("0.5"), c.mpf(0)]
    coeff = {
        "minus_inv_L": jet_neg(c, invL),
        "hZ_over_2L": _jet_scale(c, jet_mul(c, jet_mul(c, [h, c.mpf(0)], z), invL),
                                  c.mpf("0.5")),
        "minus_half": jet_neg(c, half),
        "d_over_L": jet_mul(c, d, invL),
        "two_hZ2_over_L2": jet_mul(c, [2 * h, c.mpf(0)],
                                    jet_mul(c, jet_mul(c, z, z), invL2)),
        "two_Zd_over_L2": jet_mul(c, [c.mpf(2), c.mpf(0)],
                                   jet_mul(c, jet_mul(c, z, d), invL2)),
    }
    D = jet_mul(c, coeff["minus_inv_L"], macro_I1["q"])
    D = jet_add(c, D, jet_mul(c, coeff["hZ_over_2L"], macro_I1["q_Z"]))
    D = jet_add(c, D, jet_mul(c, coeff["minus_half"], macro_I1["q_squared"]))
    D = jet_add(c, D, jet_mul(c, coeff["d_over_L"], macro_I1["U_q_Z"]))
    D = jet_add(c, D, jet_mul(c, coeff["two_hZ2_over_L2"], macro_I1["U_squared"]))
    D = jet_add(c, D, jet_mul(c, coeff["two_Zd_over_L2"], macro_I1["U_U_Z"]))
    return dict(
        normalized_by_Wmacro=True,
        Wmacro=geometry["Wmacro_record"],
        positive_geometry_K=geometry["K"],
        six_integrals={name: _rows(jet) for name, jet in macro_I1.items()},
        D_rows=_rows(D),
        coefficientsC1rows={name: _rows(jet) for name, jet in coeff.items()},
        formula=(
            "D=-I[q]/L+hZ*I[q_Z]/(2L)-I[q_squared]/2+"
            "d*I[U_q_Z]/L+2*hZ^2*I[U_squared]/L^2+"
            "2*Z*d*I[U_U_Z]/L^2"
        ),
        coefficient_definitions=dict(
            h="1-delta", d="1-Z^2", L="1-delta*Z^2",
            hZ="(1-delta)*Z", L2="L^2"
        ),
        positive_geometry_multiplied_once=True,
        finite_N_correction_added=False,
        finite_N_correction_open=True,
        ordinary_Z_orders=[0, 1],
    )


def _endpoint_inclusion(parent, endpoint_guard):
    """Consume the saved s=2 endpoint/domain receipts before macro lifting."""
    ea = parent["actual_micro_Ea_mixed_bounds"]
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    ea_domain = ea["domain"]
    b2_domain = b2["domain"]
    expected_z = ["37/100", "93/250"]
    require(ea_domain.get("s") == [0, 2] and b2_domain.get("s") == [0, 2],
            "Both saved source packets must cover the full s=0..2 micro chart")
    require(ea_domain.get("Z") == expected_z and b2_domain.get("Z") == expected_z,
            "Ea and B2 Z tubes must agree with the macro Z tube")
    require(ea_domain.get("rho") == ["4", "41/10"],
            "Ea radial endpoint/domain receipt changed")
    require(b2_domain.get("physical_y") == "log(R/Ra)=hb*s",
            "B2 physical-y endpoint receipt changed")
    native = b2["native_angular_and_axial_provider"]
    require(len(native["native_phi"]) >= ORDER + 2 and
            len(native["native_U"]) >= ORDER + 1,
            "Native endpoint provider is missing its Z4 source rows")
    return dict(
        passed=True,
        Ea_domain=ea_domain,
        b2_domain=b2_domain,
        common_Z_tube=expected_z,
        endpoint_s=2,
        endpoint_included_in_Ea_and_B2=True,
        comparison_micro_end_at_s2=True,
        comparison_endpoint_source=(
            "C6 __init__.comparison_micro_end: a.at(barphi,s,c(2)), "
            "a.at(barU,s,c(2)); barphi2/barU2 consume these rows"
        ),
        actual_macro_initial_endpoint_source=(
            "C6 __init__.actual_fields initial=self.actual_micro_end"
        ),
        actual_macro_nonreset_moment_source=(
            "C6 __init__.histories initial=self.actual_moment_end; "
            "native macro moment memory is not reset"
        ),
        endpoint_guard=endpoint_guard,
    )


def _source_contract(parent):
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    ea = parent["actual_micro_Ea_mixed_bounds"]
    native = b2["native_angular_and_axial_provider"]
    core_gauge = ea["source_AST"]["core_gauge"]
    amplitude = b2["whole_tube_amplitude"]
    require(b2["passed"] and ea["passed"] and native["passed"]
            and native["ordinary_Z0_Z4"] and native["physical_y0_y1"]
            and native["native_and_comparison_never_identified"],
            "Current native C6 source packets are required")
    require(core_gauge["passed"] and
            "gradient=a.div(a.mul(L,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))"
            in core_gauge["required_text"] and
            "Gjet=c3.C3Function(G,gradient.value,gradient.Z,gradient.ZZ)"
            in core_gauge["required_text"],
            "gaxis must be tied to the current core gradient rational jets")
    require(amplitude["passed"] and
            amplitude["exact_log_F0"] == "-logC-Lambda*G; 0<=G<=Gbar on this tube" and
            amplitude["F0_squared_ordinary_Z0_Z4_coarse_absolute_cap_one"] is True and
            amplitude["selected_logC_not_replaced_by_minimum_guard"] is True,
            "Current F0/logC amplitude receipt is required")
    numeric = parent["actual_micro_background_I1_numeric"]
    require(numeric["passed"] and numeric["actual_signed_background_I1"]
            and numeric["signed_source_rows"]["signed_function_values_are_enclosures_not_cap_centers"],
            "Saved signed whole-micro native background is required")
    macro_graph = parent["actual_micro_macro_correlated_m15"]
    require(macro_graph["passed"] and macro_graph["actual_macro_six_I1_source"]["physical_measure"][
        "R_squared_length_applied_once"], "Accepted macro target graph receipt is required")
    return dict(
        passed=True,
        field_roles=dict(
            comparison_phi_U="actual_micro_Ea_mixed_bounds.comparison_source_Z0_Z4_enclosures",
            comparison_incoming_HKM="actual_micro_Ea_mixed_bounds.comparison_source_Z0_Z4_enclosures",
            comparison_incoming_ABC="actual_micro_b_p2_mixed_bounds.comparison_complete_ABC_history_rows",
            native_pre_finite_N_U_b_M_b="actual_micro_background_I1_numeric.signed_source_rows",
            native_pre_finite_N_U_b_direct=(
                "actual_micro_b_p2_mixed_bounds.native_angular_and_axial_provider.native_U"
                " (exact-equal to saved signed U rows)"
            ),
            native_pre_finite_N_M_b_endpoint=(
                "actual_micro_background_I1_numeric.signed_source_rows.M"
                " (saved signed whole-micro endpoint source)"
            ),
            native_pre_finite_N_Q_b="new same linear Q recovery from U_b/M_b",
            full_finite_N_graph_target="parent.actual_micro_macro_correlated_m15.actual_macro_six_I1_source",
            gaxis="same current core gradient rational jets (gradient.value, gradient.Z, gradient.ZZ)",
            exact_log_F0="actual_micro_b_p2_mixed_bounds.whole_tube_amplitude.exact_log_F0",
        ),
        gaxis_from_core_gradient_rational_jets=True,
        core_gradient_source_receipt=core_gauge,
        exact_log_F0=amplitude["exact_log_F0"],
        selected_logC_not_replaced_by_minimum_guard=True,
        F0_squared_cap_receipt=True,
        native_pre_finite_N_background=True,
        signed_background_values_are_enclosures=True,
        finite_N_graph_not_numerically_evaluated=True,
        finite_N_delta_U_M_to_delta_Q_reserved=True,
        full_finite_N_enclosed=False,
        fullfiniteN_enclosed=False,
        full_pressure_global_closure=False,
        outer_RQy_graph_retained_as_future_target=True,
        finite_N_radius_not_added=True,
        accepted_numeric_builder_not_called=True,
        accepted_producer_checker_owner_norm_recurrence_quadrature_not_called=True,
    )


def build(root, parent):
    """Build current native pre-finite-N macro signed C1 enclosures."""
    root = Path(root)
    disk_parent = _load(root / PARENT)
    require(sha(root / PARENT) == PARENT_SHA256 and parent == disk_parent,
            "Current correlated parent report is required")
    require(parent.get("checks", {}).get("passed") is True,
            "Current correlated parent checks must pass")
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256, passed=True)}
    for path in (B2_REPORT, C6_SOURCE, CORE_SOURCE, CORE_CHECK, JET_SOURCE,
                 SIGNED_HELPER):
        pins[path] = _pin(root, parent, path)
    guards = _source_guards(root)
    contract = _source_contract(parent)
    endpoint_inclusion = _endpoint_inclusion(parent, guards["c6_endpoint_inclusion"])
    from mpmath.ctx_iv import MPIntervalContext
    c = MPIntervalContext()
    c.prec = 600
    geometry = _macro_geometry(c, parent)
    comparison = _comparison_histories(c, parent)
    directions = _comparison_directions(c, comparison, parent, geometry)
    drive = _full_signed_axial_drive(c, comparison, directions, parent, geometry)
    fields = _actual_macro_fields(c, comparison, directions, drive, parent, geometry)
    Q = _recover_Q(c, fields, directions)
    integrands = _six_integrands(c, fields, Q)
    I1 = _macro_I1_and_D(c, integrands, directions, geometry)
    target = parent["actual_micro_macro_correlated_m15"]["actual_macro_six_I1_source"]
    signed_rows = dict(
        comparison_phi=_rows(comparison["phi"]),
        comparison_U=_rows(comparison["U"]),
        comparison_histories={name: _rows(comparison["hulls"][name])
                              for name in HISTORY_NAMES},
        native_initial_phi=_rows(fields["initial_phi"]),
        native_initial_U_b=_rows(fields["initial_U"]),
        native_initial_M_b=_rows(fields["initial_M"]),
        native_macro_phi_b=_rows(fields["phi"]),
        native_macro_U_b=_rows(fields["U"]),
        native_macro_M_b=_rows(fields["M"]),
        native_macro_Q_b=_rows(Q),
        six_density_jets={name: _rows(jet) for name, jet in integrands.items()},
        D_over_R=_rows(directions["D_over_R"]),
        full_signed_axial_drive=_rows(drive["full_signed_axial_drive"]),
    )
    return dict(
        passed=True,
        source_family=parent["source_family"],
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_pins=pins,
        source_guards=guards,
        source_contract=contract,
        endpoint_inclusion=endpoint_inclusion,
        source_functional_identities=_source_functional_identities(),
        domain=dict(
            macro_coordinate="m",
            m=["0", "1"],
            Z=["37/100", "93/250"],
            ordinary_Z_orders=[0, 1, 2, 3],
            I1_ordinary_Z_orders=[0, 1],
            pre_finite_N_native_background=True,
        ),
        geometric_measure=geometry,
        comparison=dict(
            frozen_phi_U_from_micro2=True,
            source=comparison["source"],
            abc_source=comparison["abc_source"],
            incoming_rows={name: _rows(comparison["incoming"][name])
                           for name in HISTORY_NAMES},
            equilibrium_rows={name: _rows(comparison["equilibrium"][name])
                              for name in HISTORY_NAMES},
            history_rows={name: _rows(comparison["hulls"][name])
                          for name in HISTORY_NAMES},
            history_receipts=comparison["receipts"],
            decay_is_Z_independent=True,
            rows_are_convex_hulls=True,
            no_micro_norm_copied_as_macro_value=True,
        ),
        comparison_directions=dict(
            Z=_rows(directions["Z"]),
            d=_rows(directions["d"]),
            L=_rows(directions["L"]),
            gaxis=_rows(directions["gaxis"]),
            A=_rows(directions["A_reduced"]),
            B=_rows(directions["B_reduced"]),
            W=_rows(directions["W_reduced"]),
            D_over_R=_rows(directions["D_over_R"]),
            F0_inverse_not_materialized=True,
            F0_Z_over_F0=directions["F0_Z_over_F0"],
            Lambda_not_materialized=False,
            Lambda_materialized_from_representable_log=True,
            dimensionless_exponentials_materialized=(
                "exp(-2*ell)", "exp(logLambda)", "exp(2*logPstar)"
            ),
            Lambda_definition="exp(logLambda), logLambda=4*logPstar+1000",
        ),
        full_signed_axial_drive=dict(
            hydro=_rows(drive["hydro"]),
            pressure=_rows(drive["pressure"]),
            dressed_swirl=_rows(drive["dressed_swirl"]),
            full_drive=_rows(drive["full_signed_axial_drive"]),
            P0=_rows(drive["p0"]),
            Pstar_squared=record(drive["Pstar_squared"]),
            F0_squared_cap=_rows(drive["F0_squared_cap"]),
            R_bound=_rows(drive["R_bound"]),
            R_squared_bound=_rows(drive["R_squared_bound"]),
            F0_squared_cancellation_exact=True,
            no_inverse_F0=True,
            exact_log_F0=drive["exact_log_F0"],
            log_F0=drive["log_F0"],
            selected_logC_not_replaced_by_minimum_guard=True,
            F0_squared_coarse_cap_consumed=True,
            gaxis_source="same current core gradient rational jets",
            pressure_source_units=drive["source_units"],
        ),
        actual_macro_fields=dict(
            initial_phi=_rows(fields["initial_phi"]),
            initial_U_b=_rows(fields["initial_U"]),
            initial_M_b=_rows(fields["initial_M"]),
            initial_U_source=fields["initial_U_source"],
            initial_M_source=fields["initial_M_source"],
            logphi_initial=_rows(fields["initial_logphi"]),
            logphi_actual=_rows(fields["actual_logphi"]),
            phi_b=_rows(fields["phi"]),
            U_b=_rows(fields["U"]),
            M_b=_rows(fields["M"]),
            phi_over_comparison=_rows(fields["phi_over_comparison"]),
            Q_b=_rows(Q),
            M_is_nonreset_hull=True,
            native_M_not_comparison_M=True,
            native_initial_rows_are_saved_signed_whole_micro=True,
        ),
        new_macro_error_log_proofs=fields["error_log_proofs"],
        signed_source_rows=signed_rows,
        macro_I1=I1,
        I1=I1,
        D=I1["D_rows"],
        target_full_finite_N_graph=dict(
            parent_path=PARENT,
            target_source_handles={
                "physical_fields": target["physical_fields"],
                "source": target["source"],
                "integrals": target["integrals"],
                "physical_measure": target["physical_measure"],
            },
            retained_for_future_comparison=True,
            numerically_evaluated=False,
            finite_N_delta_U_M_to_delta_Q_reserved=True,
        ),
        limits=dict(
            finite_N_correction="open",
            finite_N_correction_added=False,
            finite_N_error_radius_available=False,
            full_finite_N_enclosed=False,
            fullfiniteN_enclosed=False,
            outer_RQy_deleted=False,
            outer_RQy_future_target_retained=True,
            complete_pressure_global_closure=False,
            global_five_moment_closure=False,
            current_scope="native pre-finite-N macro six-I1 signed ordinary C1",
        ),
        actual_native_pre_finite_N_macro_I1=True,
        signed_macro_I1_enclosures=True,
        finite_N_correction_open=True,
        normalized_by_Wmacro=True,
        Wmacro_definition="Wmacro=10000",
        accepted_ancestors_replayed=False,
    )


__all__ = ["build"]
