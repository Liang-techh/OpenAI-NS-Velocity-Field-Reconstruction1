"""Build actual second-switch/post-power pressure and cumulative m15."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073'), str(ROOT/'scripts')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_r110_m15 import OriginalR110M15
from openai_ns_reconstruction.original_macro_m15 import proxy_modes
from openai_ns_reconstruction.original_first_switch_m15 import frozen_coefficients
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME = PREFIX+'compliant_current_original_r110_m15_numeric.json'
SOURCE = PREFIX+'compliant_current_original_prescribed_R100_Rm_C6_functions.py'


def original_controls():
    tree = ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'second_switch_fields']
    assert len(nodes) == 1
    s, h, z = sy.symbols('s hbS Z', positive=True)
    phi0, U0, R = 2+z+z*z, 1+z*z, 100*sy.exp(h*(1+s))
    g = SimpleNamespace(symbol=sy.Symbol, one=sy.Integer(1), zero=sy.Integer(0), constant=sy.Rational,
                        mul=lambda *v: sy.prod(v), sub=lambda a, b: a-b)
    a = SimpleNamespace(g=g, fixed=sy.sympify, add=lambda *v: sum(v), mul=lambda *v: sy.prod(v),
        scale=lambda a, b: a*sy.sympify(b), exp=sy.exp,
        at=lambda value, variable, arg: value.subs(variable, arg),
        integral=lambda value, variable, lower, upper, **kwargs: sy.Integral(value, (variable, lower, upper)))
    sig = lambda graph, phase: sy.Function('sigma')(phase)
    scope = dict(sigma=sig)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(DIRECTORY/SOURCE), 'exec'), scope)
    result = scope['second_switch_fields'](a, dict(D_over_R=1+z/3), dict(phi=phi0, U=U0), s, h, R)
    slope = -h*R*(1+z/3)*(1-sig(g, s))/2-sy.Rational(2, 5)*sig(g, s)
    count = 0
    for expr in (result['U']-U0, result['U_y'], result['phi_y']-result['phi']*slope,
                 result['angular_phase_density']-h*slope):
        for k in (0, 1):
            assert sy.cancel(sy.diff(expr, z, k)) == 0; count += 1
    # Execute the original post-switch actual-field assignment in isolation.
    matches = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name) and node.value.func.id == 'dict'):
            keywords = {kw.arg: kw.value for kw in node.value.keywords}
            if set(keywords) == {'phi', 'U', 'phi_y', 'U_y'} and isinstance(keywords['U'], ast.Subscript):
                matches.append(node)
    assert len(matches) == 1
    scope = dict(a=a, phi=phi0, incoming=dict(actual=dict(U=U0)))
    exec(compile(ast.Module(body=matches, type_ignores=[]), str(DIRECTORY/SOURCE), 'exec'), scope)
    result = scope['actual']
    assert result['U'] == U0 and result['U_y'] == 0
    assert sy.cancel(result['phi_y']+sy.Rational(2, 5)*result['phi']) == 0
    return dict(passed=True, original_second_switch_value_Z1_identities=count,
                original_post_power_constant_U_and_power_identities=3, source=SOURCE)


def independent_primitives():
    z, t, h, eps, Y = sy.symbols('Z phase hbS epsilon Y', positive=True)
    de, axis = sy.Rational(1, 17), sy.Rational(7, 5)
    psi, mean = 1+z+z*z/3+z**3/7, sy.Rational(2, 3)-z/2+z*z/5-z**3/11
    U = 4*z+axis+eps*psi
    A, D = 1-de*z*z, 1-z*z
    L = sy.log(sy.Rational(11, 10))
    c = mp.mp.clone(); c.dps = 100
    conditions = geometries = viscous = 0
    with mp.workdps(100):
        nodes, weights = mp.gauss_quadrature(48, 'legendre')
        for chart in ('second_switch', 'post_power'):
            eta = h*(1+t) if chart == 'second_switch' else 2*h+(L-2*h)*t
            length = h if chart == 'second_switch' else L-2*h
            M = U+eps*(mean-psi)*sy.exp(-Y-eta)
            Q = (2*z*U-(1-de)*z*M-D*sy.diff(M, z))/A
            Qy = sy.diff(Q, t)/length; R = 100*sy.exp(eta)
            assert sy.cancel(sy.diff(M, t)/length-(U-M)) == 0
            assert sy.simplify(sy.diff(R*Qy, t)) == 0
            assert sy.simplify(sy.diff(R*Qy, t, z)) == 0
            flow = ((Q+(1-de)*z*sy.diff(Q, z)/2+Qy)/A+Q*(Q/2+Qy)
                    +U*(-2*z*Q+D*sy.diff(Q, z)-2*z*Qy)/A)
            args = t, z, h, eps, Y
            direct = [sy.lambdify(args, sy.diff(R*R*length*flow/2, z, k), 'mpmath') for k in (0, 1)]
            end = [sy.lambdify(args, sy.diff(R*Qy, z, k), 'mpmath') for k in (0, 1)]
            for hh in (mp.mpf('.003'), mp.mpf('.02')):
                for stop in (mp.mpf('.2'), mp.mpf('.8'), mp.mpf(1)):
                    for zz in (mp.mpf('-.3'), mp.mpf('.3')):
                        def jet(expr):
                            return J(c, [c.mpf(str(sy.N(sy.diff(expr, z, k).subs(z, sy.Rational(str(zz))), 100)))/c.factorial(k) for k in range(4)])
                        modes, _ = proxy_modes(J, c, str(zz), c.mpf(1)/17, c.mpf(7)/5, jet(psi), jet(mean))
                        H = frozen_coefficients(J, c, str(zz), c.mpf(1)/17, modes)
                        ee, yy, ll = mp.mpf('.23'), mp.mpf('1.7'), mp.log(mp.mpf(11)/10)
                        for k in (0, 1):
                            exact = mp.mpf(0)
                            for (degree, n), row in H.items():
                                rate = 2-n
                                if chart == 'second_switch':
                                    kernel = hh*stop if not rate else mp.exp(rate*hh)*mp.expm1(rate*hh*stop)/rate
                                    proxy = hh*stop
                                    bound = 0 if not rate else hh**2*rate*stop*(1+stop)*mp.exp(rate*hh*(1+stop))
                                else:
                                    kernel = (ll-2*hh)*stop if not rate else mp.exp(2*rate*hh)*mp.expm1(rate*(ll-2*hh)*stop)/rate
                                    proxy = ll*stop if not rate else mp.expm1(rate*ll*stop)/rate
                                    bound = hh*proxy*(2*rate+2/ll)*mp.exp(2*rate*hh)
                                assert abs(kernel-proxy) <= bound+mp.mpf('1e-90'); geometries += 1
                                exact += 10000*c.mpf(row[k])*ee**degree*mp.exp(-n*yy)*kernel
                            want = sum((w*direct[k]((node+1)*stop/2, zz, hh, ee, yy)*stop/2
                                        for node, w in zip(nodes, weights)), mp.mpf(0))
                            assert abs(exact-want) < mp.mpf('1e-72')*max(1, abs(want)); conditions += 1
                            lower, upper = end[k](0, zz, hh, ee, yy), end[k](stop, zz, hh, ee, yy)
                            assert lower != 0 and upper != 0
                            assert abs(upper-lower) < mp.mpf('1e-80')*max(1, abs(lower)); viscous += 1
    return dict(passed=True, independent_physical_half_S1_value_Z1_integrals=conditions,
                per_mode_native_geometry_error_conditions=geometries,
                nonzero_viscous_endpoint_cancellation_conditions=viscous,
                exact_mean_and_constant_RQy_identities_per_chart=3)


def source_mean_errors():
    c = mp.mp.clone(); c.dps = 100
    count = 0
    # Independent scalar source Psi change = h^2 B*(s-s^2/2).
    # The exact mean follows its positive native-y convolution.
    for h in (c.mpf('.003'), c.mpf('.02')):
        Lp = c.log(c.mpf(11)/10)-2*h
        B = lambda z: 1-z/3+z*z/5+z**3/9
        for q in (c.mpf('.2'), c.mpf('.8'), c.mpf(1)):
            for k in range(4):
                b = c.diff(B, c.mpf('.31'), k); P = abs(b)/2
                dU = h*h*b/2
                dM1 = b*(h*h/2-1+(h+1)*c.exp(-h))
                dM2 = dM1*c.exp(-h*q)+dU*(1-c.exp(-h*q))
                assert abs(dM2) <= h**3*(1+q)*P+c.mpf('1e-90'); count += 1
                dMsecond = dM1*c.exp(-h)+dU*(1-c.exp(-h))
                dMpost = dMsecond*c.exp(-Lp*q)+dU*(1-c.exp(-Lp*q))
                assert abs(dMpost) <= h*h*P*(2*h+1-c.exp(-Lp*q))+c.mpf('1e-90'); count += 1
    return dict(passed=True, independent_source_mean_Z0_Z3_error_conditions=count,
                source_U_error_retains_two_switch_widths=True)


def run():
    start = time.perf_counter()
    controls, independent, source_errors = original_controls(), independent_primitives(), source_mean_errors()
    print('Original constant-U source, exact mean/RQy, native primitives and source-error propagation pass', flush=True)
    second = OriginalR110M15(precision=85); post = second.post_power_field()
    audit = mp.mp.clone(); audit.dps = 190
    samples = []; row_samples = []; inlets = ordinary = serialized = zeros = endpoints = widths = seams = 0
    def check(f, Z, coordinate):
        nonlocal inlets, ordinary, serialized, zeros, endpoints, widths
        p = f.evaluate(Z, coordinate); samples.append(p)
        row = f.rows.evaluate(Z, coordinate); row_samples.append(row)
        q, _, _ = f.rows.micro._coordinate(coordinate, False)
        _, truth, _ = f.rows.background(Z, q)
        expected = 69 if f.chart == 'second_switch' else 88
        incoming = 51 if f.chart == 'second_switch' else 69
        assert f.evaluate(Z, coordinate) is p and all(p['contract'].values()) and all(row['contract'].values())
        assert len(p['m15_C0_C1_terms']) == expected and len(p['incoming_m15_terms']) == incoming
        assert p['m15_C0_C1_terms'][:incoming] == p['incoming_m15_terms']; inlets += 1
        assert p['actual_selected_integer'] == f.inlet(Z)['actual_selected_integer']
        assert row['source_cover']['generic_derivative_chart'] == f.rows.derivative_chart
        assert p['full_prefix_cover']['cumulative_history_charts'][-1] == f.chart
        for name, field in row['physical_pressure_rows'].items():
            for k, value in enumerate(field['signed_background_ordinary_Z']):
                ordinary_value = truth[name][k]*audit.factorial(k)
                assert value['source_enclosure_exact_mpf_tuples'] == [list(v) for v in ordinary_value._mpi_]; ordinary += 1
        for term in p['m15_C0_C1_terms'][incoming:]:
            for v in term['coefficient_C0_C1']:
                ends = [audit.make_mpf(tuple(t)) for t in v['source_enclosure_exact_mpf_tuples']]
                assert max(abs(e-audit.mpf(v['approximation'])) for e in ends) <= audit.mpf(v['absolute_error_upper']); serialized += 1
                if p['exact_empty_increment']:
                    assert ends[0] == ends[1] == audit.mpf(v['approximation']) == 0; zeros += 1
        for field in row['physical_pressure_rows']['U_y']['signed_background_ordinary_Z']:
            assert all(audit.make_mpf(tuple(t)) == 0 for t in field['source_enclosure_exact_mpf_tuples']); widths += 1
        names = {t['name']: t for t in p['m15_C0_C1_terms'][incoming:]}
        volume = names[f.chart+'_finite_N_volume_error']['log_scale_terms']
        if f.chart == 'second_switch':
            assert volume[0]['original_native_parameter_name'] == 'hbS'
        else:
            assert not any(a.get('source_log_atom_only') for a in volume)
        widths += 1
        if not p['exact_empty_increment']:
            lower = p['actual_lower_finite_N_endpoint_budget']
            lr = audit.make_mpf(tuple(lower['log_radius']['source_enclosure_exact_mpf_tuples'][1]))
            for name, power in (('endpoint', 2), ('viscous', 1)):
                for k in (0, 1):
                    v = lower[name][k]
                    if v['exact_zero']: continue
                    total = p['finite_N_error_logs_before_native_measure'][name][k]
                    want = audit.make_mpf(tuple(v['log_absolute_upper']['upper_exact_mpf_tuple']))+power*lr
                    got = audit.make_mpf(tuple(total['log_absolute_upper']['upper_exact_mpf_tuple']))
                    assert got >= want or abs(got-want) <= abs(got)*audit.mpf('1e-150'); endpoints += 1
            assert not row['physical_pressure_rows']['m']['finite_N_absolute_error_log_bounds'][3]['exact_zero']
        assert not p['limits']['genuine_temporal_n_ge_2']
    for Z in ('-.3', '0', '.3'):
        for f in (second, post):
            for coordinate in ('0', '.125', '.25', '.5', '.75', '1'):
                check(f, Z, coordinate)
            print('Original cumulative pressure: '+f.chart+' Z='+Z, flush=True)
        a, b = second.rows.evaluate(Z, '1'), post.rows.evaluate(Z, '0')
        for name in ('U', 'm', 'Q', 'U_y', 'Q_y'):
            assert a['physical_pressure_rows'][name]['signed_background_ordinary_Z'] == b['physical_pressure_rows'][name]['signed_background_ordinary_Z']; seams += 1
        assert post.evaluate(Z, '1')['geometry']['exact_R110_endpoint']; seams += 1
    check(second, '.3', '.875'); check(post, '.3', '.875')
    rejected = 0
    for f in (second, post):
        for args in [('1', '1'), ('nan', '.5'), ('.3', '-.01'), ('.3', '1.01'), ('.3', 'nan')]:
            try: f.evaluate(*args)
            except ValueError: rejected += 1
            else: raise AssertionError('Invalid R110 pressure query admitted')
    hashes = dict(second.hashes); hashes.update(post.hashes)
    code = dict(second.code_hashes); code.update(post.code_hashes)
    name = 'scripts/run_original_r110_m15.py'; code[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=second.family, samples=samples, pressure_row_samples=row_samples,
        input_hashes=hashes, implementation_hashes=code,
        checks=dict(passed=True, original_controls=controls, independent_native_primitives=independent,
            independent_source_error_propagation=source_errors, actual_R110_pressure_packets=len(samples),
            exact_complete_inlet_replays=inlets, ordinary_physical_pressure_rows=ordinary,
            serialized_increment_error_conditions=serialized, exact_empty_increment_rows=zeros,
            actual_lower_finite_N_endpoint_conditions=endpoints, zero_Uy_and_native_measure_conditions=widths,
            exact_source_seam_and_R110_endpoint_conditions=seams, rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start, 3),
        limits=dict(second_switch_post_power_cumulative_signed_m15_with_complete_factored_errors=True,
            resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
            independent_preRp_signed_seam=False, global_five_moment_numeric_join=False, common_global_N=False,
            exact_heat_exterior=False, stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))


def inspect(path):
    d = json.loads(path.read_text(encoding='utf-8')); assert d['checks']['passed']
    for n, h in d['input_hashes'].items(): assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest() == h, n
    for n, h in d['implementation_hashes'].items(): assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest() == h, n
    return d


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true'); parser.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args(); path = DIRECTORY/NAME
    d = inspect(path) if args.inspect_report else run()
    if args.write_report: path.write_text(json.dumps(d, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(checks=d['checks'], execution_seconds=d['execution_seconds']), indent=2), flush=True)
