"""Build original first-switch m15 with native-width pressure primitives."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073'), str(ROOT/'scripts')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_first_switch_m15 import (
    OriginalFirstSwitchM15, density_kernel, density_difference, frozen_coefficients)
from openai_ns_reconstruction.original_macro_m15 import proxy_modes, add
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from run_original_first_micro_m15 import independent_S1
from run_original_first_switch_pressure_rows import source_controls

NAME = PREFIX+'compliant_current_original_first_switch_m15_numeric.json'


def symbolic_density(z, delta, U, Q, Uy, Qy):
    A, D = 1-delta*z*z, 1-z*z
    F = Q/(2*A)+Q*Q/4-z*z*U*U/(A*A)
    G = (-Q/(2*A)+(1-delta)*z*sy.diff(Q, z)/(4*A)-Q*Q/4
         +D*U*sy.diff(Q, z)/(2*A)+(1-delta)*z*z*U*U/(A*A)
         +z*D*U*sy.diff(U, z)/(A*A))
    Fy = Qy/(2*A)+Q*Qy/2-2*z*z*U*Uy/(A*A)
    return 2*F+G+Fy


def independent_switch_integrals():
    """Direct physical density, exact mean ODE, and viscous increment."""
    z, t, alpha, eps, Y = sy.symbols('Z s alpha epsilon Y', positive=True)
    de, axis = sy.Rational(1, 17), sy.Rational(7, 5)
    psi, mean = 1+z+z*z/3+z**3/7, sy.Rational(2, 3)-z/2+z*z/5-z**3/11
    B = 1-z/3+z*z/5+z**3/9
    A, D, E = 1-de*z*z, 1-z*z, sy.exp(-alpha*t)
    u0 = 4*z+axis
    frozenU = u0+eps*psi
    frozenM = frozenU+eps*(mean-psi)*sy.exp(-Y)*E
    def pressure(U, M):
        return (2*z*U-(1-de)*z*M-D*sy.diff(M, z))/A
    frozenQ = pressure(frozenU, frozenM)
    frozenQy = sy.diff(frozenQ, t)/alpha
    flow = ((frozenQ+(1-de)*z*sy.diff(frozenQ, z)/2+frozenQy)/A
            +frozenQ*(frozenQ/2+frozenQy)
            +frozenU*(-2*z*frozenQ+D*sy.diff(frozenQ, z)-2*z*frozenQy)/A)
    # A changing source, with a nonzero initial mean. The polynomial source
    # is an independent finite-width fixture, not the original switch control.
    U = frozenU+eps*alpha**2*B*(t-t*t/2)
    M = (frozenM+eps*B*(alpha**2*(t-t*t/2)-alpha*(1-t)-1+(alpha+1)*E))
    assert sy.cancel(sy.diff(M, t)/alpha-(U-M)) == 0
    Q = pressure(U, M); Uy, Qy = sy.diff(U, t)/alpha, sy.diff(Q, t)/alpha
    K = symbolic_density(z, de, U, Q, Uy, Qy)
    sourceflow = ((Q+(1-de)*z*sy.diff(Q, z)/2+Qy)/A+Q*(Q/2+Qy)
                  +U*(-2*z*Q+D*sy.diff(Q, z)-2*z*Qy)/A)
    assert sy.cancel(K-sourceflow/2) == 0
    radius = 100*sy.exp(alpha*t)
    C = ((1-de)*z*Uy+D*sy.diff(Uy, z))/A
    identity = sy.diff(radius*Qy, t)-sy.diff(2*z*radius*Uy/A, t)+alpha*radius*C
    assert sy.cancel(identity) == 0 and sy.cancel(sy.diff(identity, z)) == 0
    args = (t, z, alpha, eps, Y)
    frozen = [sy.lambdify(args, sy.diff(10000*sy.exp(2*alpha*t)*alpha*flow/2, z, k), 'mpmath') for k in (0, 1)]
    end = [sy.lambdify(args, sy.diff(radius*Qy, z, k), 'mpmath') for k in (0, 1)]
    uyend = [sy.lambdify(args, sy.diff(2*z*radius*Uy/A, z, k), 'mpmath') for k in (0, 1)]
    integral = [sy.lambdify(args, sy.diff(alpha*radius*C, z, k), 'mpmath') for k in (0, 1)]
    kdirect = [sy.lambdify(args, sy.diff(K, z, k), 'mpmath') for k in (0, 1)]
    fields = {}
    for name, expr in dict(U=U, Q=Q, Uy=Uy, Qy=Qy).items():
        fields[name] = {d: [sy.lambdify(args, sy.diff(sy.expand(expr).coeff(eps, d), z, k)/sy.factorial(k), 'mpmath')
                            for k in range(3)] for d in (0, 1)}
    c = mp.mp.clone(); c.dps = 100
    frozen_count = viscous_count = density_count = short_count = 0
    with mp.workdps(100):
        nodes, weights = mp.gauss_quadrature(48, 'legendre')
        ee, yy = mp.mpf('.23'), mp.mpf('1.7')
        def quad(fn, zz, aa, stop):
            return sum((w*fn((node+1)*stop/2, zz, aa, ee, yy)*stop/2
                        for node, w in zip(nodes, weights)), mp.mpf(0))
        for aa in (mp.mpf('.004'), mp.mpf('.02')):
            for stop in (mp.mpf('.2'), mp.mpf('.8')):
                for zz in (mp.mpf('-.3'), mp.mpf(0), mp.mpf('.3')):
                    def jet(expr):
                        return J(c, [c.mpf(str(sy.N(sy.diff(expr, z, k).subs(z, sy.Rational(str(zz))), 100)))/c.factorial(k)
                                     for k in range(4)])
                    modes, _ = proxy_modes(J, c, str(zz), c.mpf(1)/17, c.mpf(7)/5, jet(psi), jet(mean))
                    H = frozen_coefficients(J, c, str(zz), c.mpf(1)/17, modes)
                    for k in (0, 1):
                        exact = mp.mpf(0)
                        for (d, n), h in H.items():
                            kernel = aa*stop if n == 2 else mp.expm1((2-n)*aa*stop)/(2-n)
                            exact += 10000*c.mpf(h[k])*ee**d*mp.exp(-n*yy)*kernel
                        want = quad(frozen[k], zz, aa, stop)
                        assert abs(exact-want) < mp.mpf('1e-75')*max(1, abs(want)); frozen_count += 1
                        direct = -(end[k](stop, zz, aa, ee, yy)-end[k](0, zz, aa, ee, yy))
                        reduced = -(uyend[k](stop, zz, aa, ee, yy)-uyend[k](0, zz, aa, ee, yy))+quad(integral[k], zz, aa, stop)
                        assert abs(direct-reduced) < mp.mpf('1e-75')*max(1, abs(direct)); viscous_count += 1
                    args0 = (stop, zz, aa, ee, yy)
                    rows = {name: {(d, 0): J(c, [c.mpf(fn(*args0)) for fn in fns]) for d, fns in powers.items()}
                            for name, powers in fields.items()}
                    value = density_kernel(J, c, str(zz), c.mpf(1)/17, rows)
                    for k in (0, 1):
                        got = sum(c.mpf(v[k])*ee**d for (d, n), v in value.items())
                        want = kdirect[k](*args0)
                        assert abs(got-want) < mp.mpf('1e-75')*max(1, abs(want)); density_count += 1
        # These are checks on positive finite-width remainders; no actual
        # unmaterializable original hb is replaced by a numerical fixture.
        for aa in (mp.mpf('1e-45'), mp.mpf('.004'), mp.mpf('.02')):
            for stop in (mp.mpf('.2'), mp.mpf(1)):
                for n in range(3):
                    x = aa*stop
                    exact = x if n == 2 else mp.expm1((2-n)*x)/(2-n)
                    bound = 0 if n == 2 else aa**2*stop**2*(2-n)*mp.exp((2-n)*x)/2
                    assert abs(exact-x) <= bound+mp.mpf('1e-99'); short_count += 1
    return dict(passed=True, exact_native_density_and_mean_ODE_identities=2,
        exact_viscous_value_Z1_identities=2, independent_frozen_pressure_value_Z1_integrals=frozen_count,
        independent_viscous_value_Z1_integrals=viscous_count,
        changing_source_density_value_Z1_comparisons=density_count, short_width_and_resonance_conditions=short_count)


def independent_density_errors():
    """Ordinary C1 products, including the quadratic residual-width terms."""
    z, eps, h = sy.symbols('Z epsilon h')
    signs = sy.symbols('su sq suy sqy')
    de = sy.Rational(1, 17)
    source = dict(U=1+z+eps*(sy.Rational(1, 2)+z*z+z**3/4),
                  Q=2-z/3+eps*(1+z+z*z),
                  Uy=eps*(sy.Rational(1, 3)+z/5+z*z), Qy=eps*(-sy.Rational(1, 2)+z*z/3))
    e = dict(U=(1+z+z*z)/100, Q=(2-z+z*z+z**3)/100,
             Uy=(1+z/2+z*z)/100, Qy=(3-z+z*z)/100)
    changed = {name: row+h*eps*sgn*e[name] for (name, row), sgn in zip(source.items(), signs)}
    baseK = symbolic_density(z, de, **source)
    delta = sy.expand((symbolic_density(z, de, **changed)-baseK)/h)
    fns = [sy.lambdify((z, eps, h, *signs), sy.diff(delta, z, k), 'mpmath') for k in (0, 1)]
    c = mp.mp.clone(); c.dps = 190
    iv = MPIntervalContext(); iv.dps = 130
    zz, ee, cap = c.mpf('.31'), c.mpf('.23'), iv.mpf('.01')
    def jet(expr, symmetric=False):
        values = [c.mpf(str(sy.N(sy.diff(expr, z, k).subs(z, sy.Rational(31, 100)), 190)))/c.factorial(k) for k in range(3)]
        return J(iv, [iv.mpf([-abs(v), abs(v)]) if symmetric else iv.mpf(v) for v in values])
    fields = {name: {(d, 0): jet(sy.expand(expr).coeff(eps, d)) for d in (0, 1)} for name, expr in source.items()}
    errors = {name: {(1, 0): jet(expr, True)} for name, expr in e.items()}
    result = density_difference(J, iv, '.31', iv.mpf(1)/17, fields, errors, cap)
    bound = sum((row*iv.mpf(ee)**d for (d, n), row in result.items()), J.constant(iv, 0, 1))
    broadened = {name: add(fields[name], {key: row*cap for key, row in errors[name].items()}) for name in fields}
    assert all(row.order == 2 for powers in broadened.values() for row in powers.values())
    assert all(row.order == 1 for row in density_difference(J, iv, '.31', iv.mpf(1)/17, broadened, errors, cap).values())
    count = 0
    with mp.workdps(190):
        for width in (c.mpf('.001'), c.mpf('.01')):
            for variant in itertools.product((-1, 1), repeat=4):
                for k in (0, 1):
                    want = c.mpf(fns[k](zz, ee, width, *variant))
                    lo, hi = _bounds(c, bound[k])
                    assert lo <= want <= hi; count += 1
    return dict(passed=True, independent_ordinary_C0_C1_density_error_conditions=count,
                source_baseline_keeps_Z2_before_C1_differentiation=True)


def run():
    start = time.perf_counter()
    source, controls = independent_S1(), source_controls()
    integrals, products = independent_switch_integrals(), independent_density_errors()
    print('Original pressure/control bodies, short primitives, exact viscous increment and C1 error products pass', flush=True)
    f = OriginalFirstSwitchM15(precision=85)
    audit = mp.mp.clone(); audit.dps = 190
    samples = []; inlet = zeros = serialized = atoms = histories = endpoints = terminal = 0
    def check(Z, coordinate):
        nonlocal inlet, zeros, serialized, atoms, histories, endpoints, terminal
        p = f.evaluate(Z, coordinate); samples.append(p)
        assert f.evaluate(Z, coordinate) is p
        assert all(p['contract'].values()) and len(p['m15_C0_C1_terms']) == 51
        assert p['m15_C0_C1_terms'][:28] == p['incoming_macro_m15_terms']; inlet += 1
        assert p['actual_selected_integer'] == f.inlet(Z)['actual_selected_integer']
        assert p['full_prefix_cover']['cumulative_history_charts'] == ['first_micro', 'second_micro', 'frozen_macro', 'first_switch']; histories += 1
        mean_decay = audit.make_mpf(tuple(p['inherited_core_mean_log_decay']['source_enclosure_exact_mpf_tuples'][1]))
        assert p['inherited_core_mean_log_decay'] == f.o._summary(-f.rows.micro.Ypoint, -f.rows.micro.Y)
        assert mean_decay < 0; histories += 1
        new = p['m15_C0_C1_terms'][28:]
        names = {term['name']: term for term in new}
        assert len(names) == 23
        for term in new:
            for row in term['coefficient_C0_C1']:
                ends = [audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]
                assert max(abs(x-audit.mpf(row['approximation'])) for x in ends) <= audit.mpf(row['absolute_error_upper']); serialized += 1
                if p['exact_empty_increment']:
                    assert ends[0] == ends[1] == audit.mpf(row['approximation']) == 0; zeros += 1
            for atom in term['log_scale_terms']:
                if atom.get('source_log_atom_only'):
                    assert atom['original_native_parameter_name'] in ('hbS', 'hbB')
                    assert atom['distinct_graph_nodes_with_checked_common_numeric_source']; atoms += 1
        for d in (1, 2):
            term = names[f'first_switch_inlet_error_epsilon_{d}']
            assert [atom['original_native_parameter_name'] for atom in term['log_scale_terms'][1:]] == ['hbS', 'hbB']; atoms += 1
        for d in range(3):
            for n in range(d+1):
                term = names[f'first_switch_half_S1_epsilon_{d}_mean_{n}']
                expected = f.o._summary(2*f.ctx.log(100)-d*f.o.core.logLambda-n*f.rows.micro.Ypoint,
                                       2*f.iv.log(100)-d*f.o.core.ll_iv-n*f.rows.micro.Y)
                assert term['log_scale_terms'][0] == expected; atoms += 1
                assert term['actual_function_binding']['inherited_exp_minus_nY_atom_retained']
        if not p['exact_empty_increment']:
            for name, power in (('endpoint', 2), ('viscous', 1)):
                for k in (0, 1):
                    lower = p['actual_lower_finite_N_endpoint_budget'][name][k]
                    if lower['exact_zero']: continue
                    total = p['finite_N_error_logs_before_volume_width'][name][k]
                    want = audit.make_mpf(tuple(lower['log_absolute_upper']['upper_exact_mpf_tuple']))+power*audit.log(100)
                    got = audit.make_mpf(tuple(total['log_absolute_upper']['upper_exact_mpf_tuple']))
                    assert got >= want or abs(got-want) <= abs(got)*audit.mpf('1e-150'); endpoints += 1
            term = names['first_switch_finite_N_volume_error']
            assert term['log_scale_terms'][0]['original_native_parameter_name'] == 'hbS'; atoms += 1
        if coordinate == '1':
            for row in names['first_switch_viscous_upper']['coefficient_C0_C1']:
                assert all(audit.make_mpf(tuple(t)) == 0 for t in row['source_enclosure_exact_mpf_tuples']); terminal += 1
            term = names['first_switch_finite_N_viscous_error']
            assert len(term['log_scale_terms']) == 1 and term['log_scale_terms'][0]['positive_error_bound_atom_only']; terminal += 1
        assert not p['limits']['genuine_temporal_n_ge_2'] and not p['limits']['full_corrected_NS']
    for Z in ('-.3', '0', '.3'):
        for coordinate in ('0', '.125', '.25', '.5', '.75', '1'):
            check(Z, coordinate)
        print('Actual first-switch cumulative pressure complete: Z='+Z, flush=True)
    check('.3', '.875')
    rejected = 0
    for args in [('1', '1'), ('nan', '.5'), ('.3', '-.01'), ('.3', '1.01'), ('.3', 'nan')]:
        try: f.evaluate(*args)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid first-switch moment query admitted')
    code = dict(f.code_hashes)
    for name in ('scripts/run_original_first_micro_m15.py', 'scripts/run_original_first_switch_pressure_rows.py', 'scripts/run_original_first_switch_m15.py'):
        code[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family, samples=samples, input_hashes=f.hashes, implementation_hashes=code,
        checks=dict(passed=True, original_source_S1=source, original_switch_controls=controls,
            independent_first_switch_pressure=integrals, independent_factored_density_errors=products,
            actual_first_switch_packets=len(samples), exact_macro28_inlet_replays=inlet,
            serialized_increment_error_conditions=serialized, exact_empty_switch_rows=zeros,
            native_width_and_inherited_mean_atom_conditions=atoms, complete_prefix_mean_history_conditions=histories,
            actual_lower_finite_N_pressure_endpoint_conditions=endpoints, actual_terminal_Uy_and_retained_error_conditions=terminal,
            rejected_invalid_queries=rejected), execution_seconds=round(time.perf_counter()-start, 3),
        limits=dict(first_switch_cumulative_signed_m15_with_complete_factored_errors=True,
            resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
            second_switch_post_power_m15=False, independent_preRp_signed_seam=False,
            global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
            stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))


def inspect(path):
    d = json.loads(path.read_text(encoding='utf-8')); assert d['checks']['passed']
    for n, h in d['input_hashes'].items(): assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest() == h, n
    for n, h in d['implementation_hashes'].items(): assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest() == h, n
    return d


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write-report', action='store_true'); p.add_argument('--inspect-report', action='store_true')
    args = p.parse_args(); path = DIRECTORY/NAME
    d = inspect(path) if args.inspect_report else run()
    if args.write_report: path.write_text(json.dumps(d, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(checks=d['checks'], execution_seconds=d['execution_seconds']), indent=2), flush=True)
