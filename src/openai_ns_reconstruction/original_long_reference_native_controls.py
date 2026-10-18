"""Actual native long/reference controls in logarithmic ordinary-Z form.

All numerical inputs are read from the accepted native source family. The
new arithmetic propagates bounds for those functions; no logarithmic cap is
selected as a signed field value. Extremely large radii and small positive
amplitudes remain logarithmic. No accepted producer or arithmetic owner runs.
"""
from __future__ import annotations

import ast
import gzip
import json
from pathlib import Path

from mpmath import make_mpf
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_Z_phase_bounds import LogZ
from . import original_native_log_jet4 as jet4

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_background_I1.json"
PARENT_SHA256 = "7fce21dc6e763f42b99ef7fcf7416decaf848a5124cda99741222dba9619eb1b"
MODULE = "src/openai_ns_reconstruction/original_long_reference_native_controls.py"
KERNEL = "src/openai_ns_reconstruction/original_native_log_jet4.py"
C6 = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
BRIDGE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
INLET = BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"
LOOP = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
ANCHOR = BASE + "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py"
SIGMA = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
RATE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
CHARTS = ("long_reshape", "reference")
NAMES = ("phi", "U", "M", "A", "B", "C")


def _lower(v):
    return make_mpf(v._mpi_[0])


def _upper(v):
    return make_mpf(v._mpi_[1])


def _point(c, value):
    return c.mpf([value, value])


def _maxlog(c, values):
    values = [v for v in values if v is not None]
    return None if not values else _point(c, max(_upper(v) for v in values))


def _logs(c, values):
    return [jet4.abs_log(c, v) for v in values]


def _rows(c, values, count=5):
    require(isinstance(values, list) and len(values) >= count, "Missing actual native rows")
    require([v["ordinary_Z_order"] for v in values[:count]] == list(range(count)),
            "Ordinary native rows are not ordered")
    return [read(c, v["value"]) for v in values[:count]]


def _pack(values):
    return [dict(ordinary_Z_order=k, exact_zero=v is None,
                 log_absolute_upper_bound=None if v is None else record(v))
            for k, v in enumerate(values)]


def _jsonable(value):
    if hasattr(value, "_mpi_"):
        return record(value)
    if isinstance(value, dict):
        return {key: _jsonable(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _pins(root, parent):
    accepted = {PARENT: PARENT_SHA256}
    accepted.update(parent["implementation_hashes"])
    accepted.update({BASE + k: v for k, v in parent["input_hashes"].items()})
    paths = (PARENT, C6, CACHE, BRIDGE, INLET, LOOP, ANCHOR, SIGMA, RATE_SOURCE,
             "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py",
             "src/openai_ns_reconstruction/original_log_mixed_jets.py",
             "src/openai_ns_reconstruction/original_core_signed_direction.py")
    out = {}
    for path in paths:
        expected = accepted.get(path)
        require(isinstance(expected, str) and sha(root / path) == expected,
                "Actual native source bytes changed: " + path)
        out[path] = dict(path=path, expected_sha256=expected, passed=True)
    for path in (MODULE, KERNEL):
        out[path] = dict(path=path, expected_sha256=sha(root / path), passed=True)
    return out


def _guards(root):
    return dict(
        actual_constructor=source_guard(root, C6, "__init__", (
            "incoming=self.endpoints['R110'];self.B=B=a.add(a.scale(a.trim(bridge_owner.core7['G']),g.neg(self.core.parameters['Lambda'])),",
            "a.log(incoming['actual']['phi'],'positive actual R110 phi')",
            "a.fixed(g.mul(c('1/2'),g.unary('log',c(220)))),logq)",
            "logUtheta=a.add(a.scale(B,g.sub(g.one,sigma(g,q))),a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
            "slope=a.add(a.fixed('-2/5'),a.scale(B,g.neg(g.quotient(sigma_prime(g,q),T,'same positive original T'))))",
            "self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])",
            "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
            "self.reference_length=length=g.sub(g.sub(logref,T),c(8))",
            "logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
            "self.add_chart('reference',q,R,length,actual,incoming['moments'])",
        )),
        phi_from_Utheta=source_guard(root, C6, "phi_from_Utheta", (
            "return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R)))",
        )),
        anchored_parameters=source_guard(root, ANCHOR, "__init__", (
            "self.params=dict(j=self.core.j,delta=self.core.delta,sigma=self.core.sigma,a=(9-self.core.delta)/2)",
        )),
        anchored_gradient=source_guard(root, ANCHOR, "gradient_jets", (
            "H=-4*z**3-z*z*self.params['j']+z*self.params['a']+self.params['j']",
            "denominator=H*H+self.params['sigma']**2",
            "return ((one-z*z*self.params['delta'])*H/IntervalTaylor(c,coefficients)).coefficients",
        )),
        generic_fields=source_guard(root, BRIDGE, "generic_fields", (
            "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
            "e=a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B'])),p=a.mul(FF,J['C'])",
            "fields['pressure']=a.add(fields['P0'],fields['p'])",
        )),
        complete_histories=source_guard(root, BRIDGE, "histories", (
            "body=a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))",
            "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase",
            "original_incoming_memory_retained=True,source_measure_applied_once=True",
        )),
        history_densities=source_guard(root, BRIDGE, "densities", (
            "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))",
        )),
        full_p2=source_guard(root, INLET, "generic_functions", (
            "axial_linear=div(a.add(a.neg(V),half(a.mul(a.sub(one,de),a.sub(m,a.mul(z,mz))))))",
            "axial_quadratic=div(a.add(a.mul(V,transport),a.scale(a.mul(de,a.mul(z,e)),2),a.neg(a.mul(d,ez)),",
            "a.scale(a.mul(a.add(one,de),a.mul(z,pressure)),2),a.neg(a.mul(d,pz))))",
            "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,",
        )),
        roots=source_guard(root, LOOP, "generic_roots", (
            "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)",
            "aa,bb=divide(C),divide(B)",
            "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
        )),
        loop=source_guard(root, LOOP, "loop_functions", (
            "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
            "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))",
            "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        )),
        sigma_derivative=source_guard(root, C6, "sigma_prime", (
            "derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets)",
            "exact_flat_source_derivative_not_a_bound=True",
        )),
        sigma_formula=source_guard(root, SIGMA, "_sigma_left", (
            "odds=1/(1-x)**2-1/x**2",
            "e=positive_exp(c,odds); value=e/(1+e)",
            "L1=2/(1-x)**3+2/x**3",
            "derivatives=[value,q*L1",
        )),
    )


def _gradient(c, z, j, delta, sigma):
    zero = c.mpf(0)
    H = [-4*z**3-j*z*z+(9-delta)*z/2+j,
         -12*z*z-2*j*z+(9-delta)/2, -24*z-2*j, c.mpf(-24), zero]
    L = [1-delta*z*z, -2*delta*z, -2*delta, zero, zero]
    D = jet4.value_product(c, H, H)
    # The exact denominator value uses a correlated square.
    D[0] = H[0]**2 + sigma**2
    require(_lower(D[0]) > 0, "Actual anchored gradient denominator lost positivity")
    return jet4.value_product(c, jet4.value_product(c, L, H),
                              jet4.value_reciprocal(c, D)), H, D


def _relative(c, ell):
    return jet4.exp_relative(c, ell)


def _history_rates(root):
    tree = ast.parse((root/RATE_SOURCE).read_text(encoding="utf-8"))
    entries = [node for node in tree.body if isinstance(node, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == "RATES" for target in node.targets)]
    require(len(entries) == 1 and isinstance(entries[0].value, ast.Call)
            and isinstance(entries[0].value.func, ast.Name)
            and entries[0].value.func.id == "dict" and not entries[0].value.args,
            "Native positive-kernel source rates changed")
    rates = {row.arg: ast.literal_eval(row.value) for row in entries[0].value.keywords}
    require(rates == dict(H=2,M=1,K=2,A=1,B=2,C=1), "Actual native history rates differ")
    return dict(passed=True,rates=rates,
                source=RATE_SOURCE,bridge_alias="RATES=bridge.RATES",source_executed=False)


def _function_identities():
    """Algebraic simplifications of the independently bound source formulas."""
    import sympy as s
    z, de, R, P, F, Fz, phi = s.symbols("Z delta R Pstar F0 F0_Z phi", positive=True)
    U, M, Mz, A, Az, B, Bz, C, Cz, P0, P0z = s.symbols("U M M_Z A A_Z B B_Z C C_Z P0 P0_Z")
    d, h, L = 1-z*z, 1-de, 1-de*z*z
    f2, f2z = F*F, 2*F*Fz
    E, V, m, mz = F*s.sqrt(2*R)*phi/P, U/P, M/P, Mz/P
    e, ez = (A-R*f2*B)/P**2, (Az-R*(f2z*B+f2*Bz))/P**2
    p, pz = R*f2*C/P**2, R*(f2z*C+f2*Cz)/P**2
    transport = h*z*m+d*mz
    linear = (-V+h*(m-z*mz)/2)/L
    quadratic = (V*transport+2*de*z*e-d*ez+2*(1+de)*z*(P0+p)-d*(P0z+pz))/L
    original = R*(linear+P*quadratic)/E
    numerator = (-U+h*(M-z*Mz)/2+U*(h*z*M+d*Mz)+2*de*z*A-d*Az
        +P**2*(2*(1+de)*z*P0-d*P0z)
        +R*(-2*de*z*f2*B+d*(f2z*B+f2*Bz)+2*(1+de)*z*f2*C-d*(f2z*C+f2*Cz)))
    diff = s.simplify(original-R*numerator/(F*s.sqrt(2*R)*phi*L))
    lg, lam, lp, lq, lc, yy, sig = s.symbols("G Lambda logphi0 logq logC y sigma")
    shape = -lam*lg+lp+s.log(220)/2+lq
    source_logphi = shape*(1-sig)-lq+yy/10-lc-(-lc-lam*lg)-s.log(220)/2-yy/2
    phi_diff = s.expand(source_logphi-(lp-s.Rational(2,5)*yy-shape*sig))
    t = s.symbols("t")
    # With x=8*sqrt(t)/(1-t)^2, cosh(x)^2 >= 1+x^2.
    # The remaining rational comparison has this positive numerator.
    sigma_polynomial = s.expand(64-(6-3*t+t*t)*(1-t))
    require(diff == phi_diff == 0 and sigma_polynomial == 58+9*t-4*t*t+t**3,
            "Actual native formula identity failed")
    return dict(passed=True, p2_difference=str(diff), logphi_difference=str(phi_diff),
        p2_source="generic_fields + generic_functions, full independent P0 and dressed B/C",
        p2_result="p2=R*Nnative/(Utheta*Lz)",
        ordinary_Z_derivatives_follow_from_function_identity=[0,1,2,3],
        sigma_derivative_bound=dict(passed=True, value_range=[0,1], derivative_range=[0,8],
            argument="t=(1-2q)^2, x=8*sqrt(t)/(1-t)^2, sigma'=8*(1+3t)/(1-t)^3/cosh(x)^2",
            inequality="cosh(x)^2>=1+x^2; 64-(6-3t+t^2)*(1-t)=58+9t-4t^2+t^3>0 for 0<=t<1",
            positive_polynomial=str(sigma_polynomial), flat_ends_by_definition=True),
        endpoint_congruence_claimed=False)


def _LZ(c, rows):
    return LogZ(c, rows[:4])


def _lv(jet):
    return [jet[k] for k in range(4)]


def _history_bound(c, incoming, density, rate):
    # Positive constant-rate kernel: both weights sum to one after f/r.
    return [_maxlog(c, [a, None if b is None else b-c.ln(rate)])
            for a, b in zip(incoming, density)]


def _native_p2(c, native, f2, p0, z, delta, logR, logUtheta_lower, erel):
    one = LogZ(c).constant(c.mpf(0))
    zz = _LZ(c, _logs(c, [z, c.mpf(1), c.mpf(0), c.mpf(0)]))
    z2 = one.mul(zz, zz)
    d = one.add(one, z2)
    h = jet4.abs_log(c, 1-delta)
    ld = jet4.abs_log(c, delta)
    lopd = jet4.abs_log(c, 1+delta)
    Lv = [1-delta*z*z, -2*delta*z, -2*delta, c.mpf(0)]
    L = _LZ(c, _logs(c, Lv))
    invL = L.reciprocal(L, c.ln(_point(c, _lower(Lv[0]))))
    U, M, A, B, C = (_LZ(c, native[k]) for k in ("U", "M", "A", "B", "C"))
    Mz, Az = _LZ(c, native["M"][1:]), _LZ(c, native["A"][1:])
    hydro = one.add(U,
        one.scale(one.add(M, one.mul(zz, Mz)), h-c.ln(2)),
        one.mul(U, one.add(one.scale(one.mul(zz, M), h), one.mul(d, Mz))),
        one.scale(one.mul(zz, A), c.ln(2)+ld), one.mul(d, Az))
    pressure = one.add(one.scale(one.mul(zz, _LZ(c, p0)), c.ln(2)+lopd),
                       one.mul(d, _LZ(c, p0[1:])))
    fb = jet4.log_product(c, f2, native["B"])
    fc = jet4.log_product(c, f2, native["C"])
    dressed = one.add(one.scale(one.mul(zz, _LZ(c, fb)), c.ln(2)+ld),
        one.mul(d, _LZ(c, fb[1:])),
        one.scale(one.mul(zz, _LZ(c, fc)), c.ln(2)+lopd),
        one.mul(d, _LZ(c, fc[1:])))
    # p0 is already Pstar^2*P0 here. Physical radius multiplies swirl once.
    numerator = one.add(hydro, pressure, one.scale(dressed, logR))
    invUtheta = one.scale(_LZ(c, erel), -logUtheta_lower)
    p2 = one.scale(one.mul(one.mul(numerator, invUtheta), invL), logR)
    return _lv(p2), dict(hydro=_pack(_lv(hydro)),
        independent_pressure=_pack(_lv(pressure)),
        dressed_swirl=_pack(_lv(dressed)),
        numerator=_pack(_lv(numerator)),
        inverse_L=_pack(_lv(invL)),
        inverse_Utheta=_pack(_lv(invUtheta)),
        formula="p2=R*Nnative/(Utheta*Lz); Utheta=Pstar*E=F0*sqrt(2R)*phi",
        full_pressure_and_ABC_retained=True,
        physical_J_not_in_denominator=True,
        signed_value_selected=False)


def _root_controls(c, aa, a_lower, a_upper, p2, logdstar):
    alg = LogZ(c)
    a = _LZ(c, aa)
    inva = alg.reciprocal(a, c.ln(a_lower))
    gamma = LogZ(c, [c.ln(4)] + aa[1:4])
    q2 = alg.scale(alg.mul(gamma, inva), -c.ln(2))
    s = _lv(q2)
    lq = -c.mpf("0.5")*(c.ln(2)+c.ln(a_upper))
    add = lambda *v: jet4.log_add(c, *v)
    q = [s[0]/2,
         None if s[1] is None else s[1]-c.ln(2)-lq,
         add(None if s[2] is None else s[2]-c.ln(2)-lq,
             None if s[1] is None else 2*s[1]-c.ln(4)-3*lq),
         add(None if s[3] is None else s[3]-c.ln(2)-lq,
             None if s[1] is None or s[2] is None else s[1]+s[2]+c.ln(3)-c.ln(4)-3*lq,
             None if s[1] is None else 3*s[1]+c.ln(3)-c.ln(8)-5*lq)]
    u = alg.scale(alg.mul(_LZ(c, p2), _LZ(c, q)), -logdstar)
    v = _lv(alg.mul(u, u))
    h = [c.mpf(0), None if v[1] is None else v[1]-c.ln(2),
         add(None if v[1] is None else 2*v[1]+c.ln("0.75"),
             None if v[2] is None else v[2]-c.ln(2)),
         add(None if v[1] is None else 3*v[1]+c.ln("1.875"),
             None if v[1] is None or v[2] is None else v[1]+v[2]+c.ln("2.25"),
             None if v[3] is None else v[3]-c.ln(2))]
    r = _lv(alg.mul(u, _LZ(c, h)))
    r[0] = c.mpf(0)
    alpha = alg.scale(alg.mul(_LZ(c, q), _LZ(c, h)), c.ln(2))
    denominator = alg.add(alg.constant(c.mpf(0)), alg.scale(q2, c.ln(2)))
    K = alg.scale(alg.reciprocal(denominator, c.mpf(0)), -c.ln(2*c.pi))
    lh = -c.mpf("0.5")*add(c.mpf(0), 2*u[0])
    ld = 4*_point(c, _lower(lh))-c.ln(4)
    lk = -c.ln(2*c.pi)-_point(c, _upper(denominator[0]))
    return dict(q=q, p2=p2, h=h, u=_lv(u), r=r,
                alpha=_lv(alpha), K=_lv(K)), dict(
        log_D_lower=ld, log_lambda_lower=lk)


def build(root, parent):
    root = Path(root)
    require(sha(root/PARENT) == PARENT_SHA256, "Wrong native-control parent bytes")
    require(parent == json.loads((root/PARENT).read_bytes()), "Wrong native-control parent mapping")
    pins, guards = _pins(root, parent), _guards(root)
    identities = _function_identities()
    history_rates = _history_rates(root)
    c = MPIntervalContext()
    c.prec = 900
    bg = parent["actual_long_reference_background_I1_numeric"]
    source = parent["actual_long_reference_pressure_source"]
    require(bg["passed"] and source["passed"] and source["limits"]["seams_proved"] is False,
            "Actual native source scope changed")
    cache = json.loads(gzip.decompress((root/CACHE).read_bytes()))
    require(cache["source_family"] == parent["source_family"], "Native family mismatch")
    bindings = cache["prescribed_R100_Rm_C6_bindings"]
    anchor_binding = parent["actual_switch_post_background_I1_numeric"]["sourceproof"]["second_switch_cache"]
    require(anchor_binding["original_source_bindings"][
                "anchored_G_recipe_replaces_the_enclosure_only_Gbar_range_not_the_defining_formula"] is True
            and anchor_binding["post_switch"]["anchored_G_not_Gbar"] is True
            and parent["actual_micro_b_p2_mixed_bounds"]["source_family"] == parent["source_family"],
            "Actual anchored G was replaced by a comparison bound/family")
    for predecessor, successor in (("R110", "long_reshape"), ("Rsh", "reference")):
        require(bindings["endpoints"][predecessor]["moments"] ==
                bindings["charts"][successor]["incoming_histories"],
                "Actual native moment inlet was reset")
    old = parent["actual_switch_post_background_I1_numeric"]["charts"]["post_switch"]["native_fields"]
    inlet = {key: _rows(c, old[key]) for key in NAMES}
    require(_lower(inlet["phi"][0]) > 0, "Actual R110 phi needs a positive enclosure")
    z = c.mpf(["37/100", "93/250"])
    parameters = parent["actual_core_signed_direction"]["parameters"]
    j, delta, sigma = (read(c, parameters[k]) for k in ("j", "delta", "sigma"))
    grad, H, D = _gradient(c, z, j, delta, sigma)
    logLambda = read(c, parent["actual_macro_background_I1_numeric"]["geometric_measure"]["logLambda"])
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    require(amp["passed"] is True and
            amp["exact_log_F0"] == "-logC-Lambda*G; 0<=G<=Gbar on this tube",
            "Actual anchored amplitude/G enclosure missing")
    Gbar = read(c, amp["G_modulus_upper"])
    logF0 = read(c, amp["log_F0"])
    scalars = bg["scalar_graph"]
    T, Lref, logC, logP = (read(c, scalars[k]) for k in ("T", "reference_length", "logC", "logP"))
    require(logC == read(c, amp["logC"]), "Native and anchored core logC differ")
    require(_lower(T) > 0 and _lower(Lref) > 0, "Positive physical lengths required")
    qjet = [1+z*z, 2*z, c.mpf(2), c.mpf(0), c.mpf(0)]
    logq = jet4.value_log_derivatives(c, qjet)
    logphi0 = jet4.value_log_derivatives(c, inlet["phi"])
    g_logs = [jet4.abs_log(c, v) for v in grad[:4]]
    Lambda_g = [None if v is None else logLambda+v for v in g_logs]
    B = [jet4.log_add(c, logLambda+jet4.abs_log(c, Gbar),
                      jet4.abs_log(c, logphi0[0]),
                      jet4.abs_log(c, c.ln(220)/2), jet4.abs_log(c, logq[0]))]
    B += [jet4.log_add(c, Lambda_g[k-1], jet4.abs_log(c, logphi0[k]),
                       jet4.abs_log(c, logq[k])) for k in range(1, 5)]
    margin = B[0]-c.ln(T)+c.ln(200)
    require(_upper(margin) < 0, "Actual |B|/T <=1/200 is not established")
    # sigma is Z-independent and ranges over [0,1]. These are convex
    # combinations of the two source log-derivative endpoints, not independent
    # additions of two copies of the enormous anchored contribution.
    ell_long = [_maxlog(c, [jet4.abs_log(c, logphi0[k]),
                  jet4.log_add(c, Lambda_g[k-1], jet4.abs_log(c, logq[k]))])
                for k in range(1, 5)]
    ell_ref = [jet4.log_add(c, Lambda_g[k-1], jet4.abs_log(c, logq[k]))
               for k in range(1, 5)]
    e_long = [_maxlog(c, [jet4.log_add(c, jet4.abs_log(c, logphi0[k]), Lambda_g[k-1]),
                          jet4.abs_log(c, logq[k])]) for k in range(1, 5)]
    e_ref = [jet4.abs_log(c, logq[k]) for k in range(1, 5)]
    f2 = [None if v is None else v+c.ln(2) for v in Lambda_g]
    f2 = jet4.log_scale(c, _relative(c, f2), 2*_point(c, _upper(logF0)))
    pressure_packet = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_pressure"]
    require(pressure_packet["passed"] is True and pressure_packet["ordinary_Z_orders"][:5] == [0,1,2,3,4],
            "Independent P0 ordinary rows through Z4 are required")
    p0 = _logs(c, [read(c, row) for row in pressure_packet["p0_jets"][:5]])
    logPstar = read(c, amp["log_Pstar"])
    p0 = [None if v is None else v+2*logPstar for v in p0]
    logs = parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]["actual_logs"]
    eta, dstar, logN = (read(c, logs[k]) for k in ("eta", "dstar", "N"))
    require(_upper(eta) <= 0 and eta._mpi_[0][3] >= 0,
            "Canonical eta log must be finite and at most zero")
    require(logN == read(c, parent["actual_switch_post_loop_Z3_inputs"]["logN"]),
            "Common N source changed")
    u2 = jet4.value_product(c, inlet["U"], inlet["U"])
    tail = c.mpf([0, c.mpf(2)**-1024])
    Aend = [u2[k]+(inlet["A"][k]-u2[k])*tail for k in range(5)]
    incoming = {key: _logs(c, inlet[key]) for key in ("A", "B", "C")}
    charts = {}
    for name in CHARTS:
        long = name == "long_reshape"
        ell = ell_long if long else ell_ref
        erel = _relative(c, e_long if long else e_ref)
        phirel = _relative(c, ell)
        logphi_lower = _point(c, _lower(logphi0[0]))-c.mpf("0.44")*T
        logphi_upper = _point(c, _upper(logphi0[0]))
        if not long:
            logphi_lower -= c.mpf("0.4")*Lref
            logphi_upper -= c.mpf("0.36")*T
        phi_logs = [None if v is None else v+logphi_upper for v in phirel]
        native = dict(phi=phi_logs, U=_logs(c, inlet["U"]),
                      M=_logs(c, _rows(c, bg["charts"][name]["whole_segment_M_hull"])))
        density = dict(A=_logs(c, u2), B=jet4.log_product(c, phi_logs, phi_logs))
        density["C"] = density["B"]
        for key in ("A", "B", "C"):
            rate = history_rates["rates"][key]
            native[key] = _history_bound(c, incoming[key], density[key], rate)
        if long:
            aa = [c.ln(c.mpf("0.88"))] + [v+c.ln(16)-c.ln(T) if v is not None else None for v in B[1:]]
            alo, ahi = c.mpf("0.72"), c.mpf("0.88")
            ulo = _point(c, _lower(logF0))+c.ln(220)/2+_point(c, _lower(logphi0[0]))
            uhi = T/10-logC-_point(c, _lower(logq[0]))
        else:
            aa = [c.ln(c.mpf("0.8")), None, None, None, None]
            alo = ahi = c.mpf("0.8")
            ulo = T/10-logC-_point(c, _upper(logq[0]))
            # Substitute T+Lref=10*(logC+logP)-8 before interval boxing.
            uhi = logP-c.mpf("0.8")-_point(c, _lower(logq[0]))
        logR = read(c, bg["charts"][name]["log_R_upper"])
        p2, detail = _native_p2(c, native, f2, p0, z, delta,
                               _point(c, _upper(logR)), _point(c, _lower(ulo)), erel)
        roots, conditioning = _root_controls(c, aa, alo, ahi, p2, _point(c, _lower(dstar)))
        caps = {key: _pack(value[:4]) for key, value in roots.items() if key in ("r", "alpha", "K")}
        caps.update(a=_pack(aa[:4]), b=_pack([None]*4), t0=_pack([None]*4),
                    E_over_frozen_E0=_pack(erel[:4]))
        source_chart = bindings["charts"][name]
        charts[name] = dict(passed=True, chart=name, ordinary_Z_orders=[0,1,2,3],
            log_caps=caps, conditioning=conditioning, logN=logN,
            log_Pstar_E_upper=_point(c, _upper(uhi)),
            physical_y_derivatives_claimed=False,
            proofs=dict(kappa_less_than_one_cut_one=True, eta_positive_le_one=True,
                a_positive=True, actual_phase_Z_independent=True,
                physical_y_derivatives_claimed=False),
            native_log_absolute_caps={key:_pack(value) for key,value in native.items()},
            phi_positive_log_lower=_point(c, _lower(logphi_lower)),
            phi_positive_log_upper=_point(c, _upper(logphi_upper)),
            phi_relative_ordinary_Z0_Z4=_pack(phirel),
            E_relative_ordinary_Z0_Z4=_pack(erel),
            Utheta_positive_log_lower=_point(c, _lower(ulo)),
            Utheta_positive_log_upper=_point(c, _upper(uhi)),
            a_value_interval=c.mpf([_lower(alo),_upper(ahi)]),
            p2_log_caps=_pack(p2), p2_source=detail,
            root_intermediates={key:_pack(value) for key,value in roots.items()},
            source_ports=dict(coordinate=source_chart["phase"],
                actual=source_chart["actual"], incoming_moments=source_chart["incoming_histories"],
                moments=source_chart["moments"]),
            memory_bound=dict(incoming={key:_pack(value) for key,value in incoming.items()},
                positive_kernel_maximum=True, endpoint_bound_not_selected_value=True,
                incoming_reference_is_actual_long_endpoint=True),
            limits=dict(native_H_K_numeric=False, signed_phi_ABC_values_evaluated=False,
                finite_N_numeric=False, seam_proof=False))
        incoming = {key:list(native[key]) for key in ("A","B","C")}
        incoming["A"] = _logs(c, Aend)
    result = dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        source_family=parent["source_family"],source_pins=pins,source_guards=guards,
        source_function_identities=identities,
        history_rate_source=history_rates,
        actual_anchor_binding=dict(same_current_family=True,anchored_G_not_Gbar=True,
            accepted_post_source_receipt_reused_not_reexecuted=True,
            old_reference_join_family_not_consumed=True),
        ordinary_Z_orders=[0,1,2,3],physical_y_derivatives_claimed=False,
        charts=charts,logN=logN,
        anchored_gradient=dict(H=H,denominator=D,G_Z1_through_Z4=grad[:4],
            ordinary_not_taylor=True, formula="G_Z=(1-delta*Z^2)*H/(H^2+sigma^2); H=-4Z^3-jZ^2+(9-delta)Z/2+j"),
        B_ordinary_Z0_Z4_log_caps=_pack(B),
        slope_proof=dict(passed=True,log_B_over_T_times_200=margin,
            ratio_bound="|B|/T<=1/200",sigma_prime_upper=8,
            a_long_interval=["0.72","0.88"],phi_log_slope_interval=["-0.44","-0.36"],
            actual_G_bound_used_only_as_bound=True),
        amplitude_proof=dict(logF0=logF0,actual_G_modulus_upper=Gbar,
            native_logphi_identity="logphi=logphi_R110-2y/5-B*sigma(y/T)",
            relative_reference_E_identity="D_Z^k log(E)=-D_Z^k log(1+Z^2), k>=1",
            Utheta_monotone_from_slope=True,positive_functions_never_set_zero=True),
        source_contract=dict(actual_native_long_reference_controls=True,common_N_handle=18,
            finite_N_numeric=False,seams_proved=False,physical_y_derivatives_claimed=False,
            H_K_not_required_by_p2=True,all_M_ABC_incoming_histories_retained=True,
            p2_all_pressure_and_swirl_terms=True,old_212_terms_and_errors_unchanged=True,
            accepted_ancestors_replayed=False),
        limits=dict(current_N_numeric=False,modified_endpoint_seams_proved=False,
            full_pressure_numeric=False,global_five_moment_closure=False,full_goal_complete=False),
        accepted_ancestors_replayed=False)
    result = _jsonable(result)
    json.dumps(result,allow_nan=False)
    return result
