"""Build finite original long pressure with formal signed exponential kernels."""
import argparse
import ast
import hashlib
import json
import math
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
from openai_ns_reconstruction.original_long_m15 import OriginalLongM15
from openai_ns_reconstruction.original_macro_m15 import proxy_modes
from openai_ns_reconstruction.original_first_switch_m15 import frozen_coefficients
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME = PREFIX+'compliant_current_original_long_m15_numeric.json'
SOURCE = PREFIX+'compliant_current_original_prescribed_R100_Rm_C6_functions.py'


def original_source():
    tree = ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'CurrentOriginalPrescribedR100RmC6Functions')
    init = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == '__init__')
    start = next(i for i, n in enumerate(init.body) if isinstance(n, ast.Assign)
        and isinstance(n.value, ast.Call) and any(isinstance(a, ast.Constant) and a.value == 'true_prescribed_long_reshape_fraction' for a in n.value.args))
    z, T, phi, slope = sy.symbols('Z T phi slope')
    g = SimpleNamespace(symbol=sy.Symbol, mul=lambda *v: sy.prod(v), unary=lambda name, value: sy.exp(value))
    scope = dict(g=g, c=sy.Rational, T=T)
    exec(compile(ast.Module(body=init.body[start:start+3], type_ignores=[]), str(DIRECTORY/SOURCE), 'exec'), scope)
    assert scope['y'] == T*scope['q'] and scope['R'] == 110*sy.exp(T*scope['q'])
    node = next(n for n in init.body[start:] if isinstance(n, ast.Assign)
        and any(isinstance(a, ast.Name) and a.id == 'actual' for a in n.targets))
    U = 4*z+1+z*z
    scope.update(a=SimpleNamespace(fixed=sy.sympify, mul=lambda *v: sy.prod(v)),
                 phi=phi, slope=slope, incoming=dict(actual=dict(U=U)), Utheta=sy.Symbol('ut'), logUtheta=sy.Symbol('logut'))
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(DIRECTORY/SOURCE), 'exec'), scope)
    actual = scope['actual']; conditions = 0
    for expr in (actual['U']-U, actual['U_y']):
        for k in (0, 1):
            assert sy.cancel(sy.diff(expr, z, k)) == 0; conditions += 1
    return dict(passed=True, original_long_constant_U_value_Z1_identities=conditions,
                actual_finite_T_geometry_identities=2, source=SOURCE)


def independent_primitives():
    z, t, eps, Y = sy.symbols('Z y epsilon Y')
    de, axis = sy.Rational(1, 17), sy.Rational(7, 5)
    psi, mean = 1+z+z*z/3+z**3/7, sy.Rational(2, 3)-z/2+z*z/5-z**3/11
    U = 4*z+axis+eps*psi; L = sy.log(sy.Rational(11, 10))
    M = U+eps*(mean-psi)*sy.exp(-Y-L-t)
    A, D = 1-de*z*z, 1-z*z
    Q = (2*z*U-(1-de)*z*M-D*sy.diff(M, z))/A
    Qy = sy.diff(Q, t); R = 110*sy.exp(t)
    assert sy.cancel(sy.diff(M, t)-(U-M)) == 0
    assert sy.simplify(sy.diff(R*Qy, t)) == 0
    flow = ((Q+(1-de)*z*sy.diff(Q, z)/2+Qy)/A+Q*(Q/2+Qy)
            +U*(-2*z*Q+D*sy.diff(Q, z)-2*z*Qy)/A)
    direct = [sy.lambdify((t, z, eps, Y), sy.diff(R*R*flow/2, z, k), 'mpmath') for k in (0, 1)]
    ends = [sy.lambdify((t, z, eps, Y), sy.diff(R*Qy, z, k), 'mpmath') for k in (0, 1)]
    c = mp.mp.clone(); c.dps = 100
    integrals = pairs = viscous = mass_checks = 0
    with mp.workdps(100):
        nodes, weights = mp.gauss_quadrature(48, 'legendre')
        for yy in (mp.mpf('.01'), mp.mpf('.4'), mp.mpf(2)):
            mass = -mp.expm1(-2*yy)/2
            assert mass <= min(yy, mp.mpf('.5')); mass_checks += 1
            for zz in (mp.mpf('-.3'), mp.mpf(0), mp.mpf('.3')):
                def jet(expr):
                    return J(c, [c.mpf(str(sy.N(sy.diff(expr, z, k).subs(z, sy.Rational(str(zz))), 100)))/c.factorial(k) for k in range(4)])
                modes, _ = proxy_modes(J, c, str(zz), c.mpf(1)/17, c.mpf(7)/5, jet(psi), jet(mean))
                H = frozen_coefficients(J, c, str(zz), c.mpf(1)/17, modes)
                ee, inherited = mp.mpf('.23'), mp.mpf('1.7')
                for k in (0, 1):
                    exact = mp.mpf(0)
                    for (d, n), row in H.items():
                        rate = 2-n
                        kernel = yy if not rate else mp.expm1(rate*yy)/rate
                        if rate:
                            # Only this finite independent fixture evaluates
                            # the pair. Original T/radius are never exponentiated.
                            pair = mp.exp(rate*yy)/rate-mp.mpf(1)/rate
                            assert abs(pair-kernel) < mp.mpf('1e-95'); pairs += 1
                        exact += 12100*c.mpf(row[k])*ee**d*mp.exp(-n*(inherited+mp.log(mp.mpf(11)/10)))*kernel
                    want = sum((w*direct[k]((node+1)*yy/2, zz, ee, inherited)*yy/2
                                for node, w in zip(nodes, weights)), mp.mpf(0))
                    assert abs(exact-want) < mp.mpf('1e-70')*max(1, abs(want)); integrals += 1
                    lower, upper = ends[k](0, zz, ee, inherited), ends[k](yy, zz, ee, inherited)
                    assert lower != 0 and upper != 0
                    assert abs(upper-lower) < mp.mpf('1e-85')*max(1, abs(lower)); viscous += 1
        # Finite source-error propagation from the R110 uncertainty, with
        # nonzero independent incoming U and M errors and ordinary Z0..3.
        source = 0
        P = lambda x: 1+x+x*x/3+x**3/5
        for yy in (mp.mpf('.1'), mp.mpf(2), mp.mpf(100)):
            for k in range(4):
                p = abs(mp.diff(P, mp.mpf('.31'), k))
                incoming, target = p*mp.mpf('1.04'), p
                actual = incoming*mp.exp(-yy)+target*(1-mp.exp(-yy))
                assert actual <= p*mp.mpf('2.04'); source += 1
    return dict(passed=True, independent_physical_half_S1_value_Z1_integrals=integrals,
                finite_signed_pair_kernel_comparisons=pairs, nonzero_RQy_endpoint_cancellation_conditions=viscous,
                finite_native_positive_mass_conditions=mass_checks, source_mean_Z0_Z3_propagation_conditions=source,
                exact_mean_ODE_and_RQy_identities=2)


def run():
    start = time.perf_counter(); source, independent = original_source(), independent_primitives()
    print('Original finite-T source, independent long pressure primitives and correlated kernels pass', flush=True)
    f = OriginalLongM15(precision=85); audit = mp.mp.clone(); audit.dps = 190
    samples = []; pressure = []; inlet = ordinary = serialized = zeros = pairs = resonances = endpoints = seams = memory = native = 0
    def check(Z, coordinate):
        nonlocal inlet, ordinary, serialized, zeros, pairs, resonances, endpoints, memory, native
        p = f.evaluate(Z, coordinate); samples.append(p)
        row = f.rows.evaluate(Z, coordinate); pressure.append(row)
        q, _, _ = f.rows._coordinate(coordinate, False)
        _, truth, geo = f.rows.background(Z, q)
        assert f.evaluate(Z, coordinate) is p
        assert all(p['contract'].values()) and all(row['contract'].values())
        assert len(p['m15_C0_C1_terms']) == 106 and p['m15_C0_C1_terms'][:88] == p['incoming_m15_terms']; inlet += 1
        assert p['actual_selected_integer'] == f.inlet(Z)['actual_selected_integer']
        assert row['source_cover']['generic_derivative_chart'] == 'reshape'
        assert p['geometry']['original_T_singleton_exact_mpf_tuple'] == f.rows.long.T_tuple
        assert not p['geometry']['radius_and_enormous_exponentials_materialized']
        assert not row['finite_N_axial_radial_error_components']['extra_hb_factor']; native += 1
        for name, field in row['physical_pressure_rows'].items():
            for k, value in enumerate(field['signed_background_ordinary_Z']):
                exact = truth[name][k]*math.factorial(k)
                assert value['source_enclosure_exact_mpf_tuples'] == [list(v) for v in exact._mpi_]; ordinary += 1
        for term in p['m15_C0_C1_terms'][88:]:
            for v in term['coefficient_C0_C1']:
                ends = [audit.make_mpf(tuple(t)) for t in v['source_enclosure_exact_mpf_tuples']]
                assert max(abs(e-audit.mpf(v['approximation'])) for e in ends) <= audit.mpf(v['absolute_error_upper']); serialized += 1
                if p['exact_empty_increment']:
                    assert ends[0] == ends[1] == audit.mpf(v['approximation']) == 0; zeros += 1
        names = {t['name']: t for t in p['m15_C0_C1_terms'][88:]}
        for d in range(3):
            for n in range(d+1):
                if n == 2:
                    term = names[f'long_half_S1_epsilon_{d}_mean_{n}_resonance']
                    if not p['exact_empty_increment']:
                        assert term['log_scale_terms'][1]['positive_native_log_measure_only']
                        assert term['log_scale_terms'][1]['huge_y_not_multiplied_into_coefficient']; resonances += 1
                    continue
                up = names[f'long_half_S1_epsilon_{d}_mean_{n}_upper']
                low = names[f'long_half_S1_epsilon_{d}_mean_{n}_lower']
                assert up['actual_function_binding']['correlated_coefficient_id'] == low['actual_function_binding']['correlated_coefficient_id']
                for a, b in zip(up['coefficient_C0_C1'], low['coefficient_C0_C1']):
                    aa = [audit.make_mpf(tuple(t)) for t in a['source_enclosure_exact_mpf_tuples']]
                    bb = [audit.make_mpf(tuple(t)) for t in b['source_enclosure_exact_mpf_tuples']]
                    assert aa == [-bb[1], -bb[0]]; pairs += 1
        for proof in p['geometry']['mean_memory_proofs']:
            if 'original_log_upper' in proof:
                cap = proof['original_log_upper']['log_absolute_upper']
                upper = audit.make_mpf(tuple(cap['upper_exact_mpf_tuple']))
                assert upper <= audit.log(audit.mpf('1e-1000'))
                assert proof['cap_is_bound_not_source_value']
            else:
                upper = audit.make_mpf(tuple(proof['log_magnitude_upper']['upper_exact_mpf_tuple']))
                assert upper <= -1000 and proof['exact_source_not_replaced']
            memory += 1
        if coordinate.startswith('y:'):
            assert p['coordinate']['kind'] == 'exact_decimal_R110_log_offset'
            assert geo['y']._mpi_ == f.iv.mpf(coordinate[2:])._mpi_; native += 1
        if not p['exact_empty_increment']:
            lower = p['actual_lower_finite_N_endpoint_budget']
            for name, power in (('endpoint', 2), ('viscous', 1)):
                for k in (0, 1):
                    v = lower[name][k]
                    if v['exact_zero']: continue
                    total = p['finite_N_error_logs_in_physical_units'][name][k]
                    want = audit.make_mpf(tuple(v['log_absolute_upper']['upper_exact_mpf_tuple']))+power*audit.log(110)
                    got = audit.make_mpf(tuple(total['log_absolute_upper']['upper_exact_mpf_tuple']))
                    assert got >= want or abs(got-want) <= abs(got)*audit.mpf('1e-150'); endpoints += 1
            assert not row['physical_pressure_rows']['m']['finite_N_absolute_error_log_bounds'][3]['exact_zero']
        assert not p['limits']['genuine_temporal_n_ge_2']
    for Z in ('-.3', '0', '.3'):
        for coordinate in ('0', 'y:.1', 'y:1', '.25', '.5', '1'):
            check(Z, coordinate)
        a, b = f.rows.inlet_rows.evaluate(Z, '1'), f.rows.evaluate(Z, '0')
        for name in ('U', 'm', 'Q', 'U_y', 'Q_y'):
            assert a['physical_pressure_rows'][name]['signed_background_ordinary_Z'] == b['physical_pressure_rows'][name]['signed_background_ordinary_Z']; seams += 1
        print('Original full finite long pressure and106-term m15: Z='+Z, flush=True)
    check('.3', '.875')
    rejected = 0
    for args in [('1', '.5'), ('nan', '.5'), ('.3', '-.1'), ('.3', '1.1'), ('.3', 'y:-1'), ('.3', 'y:nan')]:
        try: f.evaluate(*args)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid original long pressure query admitted')
    code = dict(f.code_hashes); name = 'scripts/run_original_long_m15.py'
    code[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family, samples=samples, pressure_row_samples=pressure,
        input_hashes=f.hashes, implementation_hashes=code,
        checks=dict(passed=True, actual_original_long_source=source, independent_long_primitives=independent,
            actual_long_packets=len(samples), exact_R110_88term_inlet_replays=inlet,
            ordinary_physical_pressure_rows=ordinary, serialized_increment_error_conditions=serialized,
            exact_empty_increment_rows=zeros, correlated_pair_coefficient_conditions=pairs,
            positive_resonant_measure_atoms=resonances, retained_lower_finite_N_endpoint_conditions=endpoints,
            exact_R110_source_seams=seams, coefficient_before_decay_cap_conditions=memory,
            native_length_and_log_offset_conditions=native, rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start, 3),
        limits=dict(long_cumulative_signed_m15_with_complete_factored_errors=True,
            resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
            later_long_Rm_pressure=False, independent_preRp_signed_seam=False,
            global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
            stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))


def inspect(path):
    data = json.loads(path.read_text(encoding='utf-8')); assert data['checks']['passed']
    for name, digest in data['input_hashes'].items(): assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest, name
    for name, digest in data['implementation_hashes'].items(): assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true'); parser.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args(); path = DIRECTORY/NAME
    data = inspect(path) if args.inspect_report else run()
    if args.write_report: path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(checks=data['checks'], execution_seconds=data['execution_seconds']), indent=2), flush=True)
