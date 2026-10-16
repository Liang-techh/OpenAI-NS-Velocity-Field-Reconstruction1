"""Build cumulative original macro m15 through physical R=100."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073'), str(ROOT/'scripts')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_macro_m15 import OriginalMacroM15, proxy_modes, proxy_integral, kernels, difference_products
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from run_original_first_micro_m15 import independent_S1

NAME = PREFIX+'compliant_current_original_macro_m15_numeric.json'


def independent_exponential_integrals():
    """Compare the analytical modes to the original physical pressure density."""
    c = mp.mp.clone(); c.dps = 90
    v, z, Y, eps = sy.symbols('v Z Y epsilon', positive=True)
    delta, axis_j = sy.Rational(1, 17), sy.Rational(7, 5)
    psi = 1+z+z*z/3+z**3/7
    mean = sy.Rational(2, 3)-z/2+z*z/5-z**3/11
    U = 4*z+axis_j+eps*psi
    M = U+eps*(mean-psi)*sy.exp(-Y*v)
    A, D = 1-delta*z*z, 1-z*z
    Q = (2*z*U-(1-delta)*z*M-D*sy.diff(M, z))/A
    Qy = sy.diff(Q, v)/Y
    R = 4*eps*sy.exp(Y*v)
    # Both finite proxy endpoints are nonzero, but their RQy product is
    # exactly constant. This does not zero either actual endpoint error.
    assert sy.cancel(sy.diff(R*Qy, v)) == 0
    assert sy.cancel(sy.diff(R*Qy, v, z)) == 0
    flow = (Q+(1-delta)*z*sy.diff(Q, z)/2+Qy)/A+Q*(Q/2+Qy)
    flow += U*(-2*z*Q+D*sy.diff(Q, z)-2*z*Qy)/A
    direct = [sy.lambdify((v, z, Y, eps), expr, 'mpmath') for expr in
              (R*R*Y*flow/2, sy.diff(R*R*Y*flow/2, z))]
    # Lambdify mpmath functions use its module context; sufficient precision
    # is held for the independent native Gauss integrations only.
    conditions = modal = 0
    with mp.workdps(90):
        nodes, weights = mp.gauss_quadrature(64, 'legendre')
        for yy in ('.3', '1.7'):
            for q in ('.2', '.8'):
                y, stop = c.mpf(yy), c.mpf(q)
                E, W = kernels(c, y, stop)
                upper = 4*c.mpf('.23')*c.exp(y*stop)
                for mode in range(3):
                    integral = sum((c.mpf(w)*y*c.exp(-(2*y)*(stop-(c.mpf(node)+1)*stop/2))
                                    *c.exp(-mode*y*(c.mpf(node)+1)*stop/2)*stop/2
                                    for node, w in zip(nodes, weights)), c.mpf(0))
                    assert abs(integral-W[mode]) < c.mpf('1e-75'); modal += 1
                for zz in ('-.3', '0', '.3'):
                    def jet(expr):
                        return J(c, [c.mpf(str(sy.N(sy.diff(expr, z, k).subs(z, sy.Rational(zz)), 90)))/c.factorial(k)
                                     for k in range(4)])
                    modes, qa = proxy_modes(J, c, zz, c.mpf(1)/17, c.mpf(7)/5, jet(psi), jet(mean))
                    parts, _, _ = proxy_integral(J, c, zz, c.mpf(1)/17, modes, y, stop)
                    for k in (0, 1):
                        got = upper*upper*sum(parts[d][k]*c.mpf('.23')**d for d in range(3))
                        want = sum((c.mpf(w)*c.mpf(direct[k]((mp.mpf(node)+1)*mp.mpf(q)/2,
                            mp.mpf(zz), mp.mpf(yy), mp.mpf('.23')))*stop/2 for node, w in zip(nodes, weights)), c.mpf(0))
                        assert abs(got-want) < c.mpf('1e-70')*max(c.mpf(1), abs(want)); conditions += 1
    # Small-prefix behavior and the exact resonant kernel do not divide by0.
    E, W = kernels(c, c.mpf('1e18'), c.mpf('1e-100'))
    assert abs(W[0]/c.mpf('1e-82')-1) < c.mpf('1e-75')
    assert W[2] > 0
    return dict(passed=True, independent_pressure_density_value_Z1_integrals=conditions,
                independent_positive_modal_integrals=modal, exact_constant_proxy_R_Qy_identities=2,
                very_small_prefix_and_resonance_conditions=2)


def independent_error_products():
    # Evaluate the independent point fixture above interval precision, so
    # equality at a directed bound is not spoiled by point rounding.
    c = mp.mp.clone(); c.dps = 190
    iv = MPIntervalContext(); iv.dps = 130
    a, b = [c.mpf('1.3'), c.mpf('-.2')], [c.mpf('.7'), c.mpf('.4')]
    ea, eb = [c.mpf('.03'), c.mpf('.04')], [c.mpf('.02'), c.mpf('.01')]
    scalar = lambda values: J(iv, [iv.mpf(v) for v in values])
    error = lambda values: J(iv, [iv.mpf([-v, v]) for v in values])
    diff = difference_products({(0, 0): scalar(a)}, {(0, 0): scalar(b)},
        {(1, 0): error(ea)}, {(1, 0): error(eb)}, iv.mpf('.01'))
    bound = diff[(1, 0)]+diff[(2, 0)]*iv.mpf('.23')
    count = geometry = 0
    for h in (c.mpf('.001'), c.mpf('.01')):
        for signs in itertools.product((-1, 1), repeat=4):
            x, y = [ea[k]*signs[k] for k in (0, 1)], [eb[k]*signs[k+2] for k in (0, 1)]
            ee = c.mpf('.23')
            got = [a[0]*y[0]+b[0]*x[0]+h*ee*x[0]*y[0],
                   a[1]*y[0]+a[0]*y[1]+b[1]*x[0]+b[0]*x[1]+h*ee*(x[1]*y[0]+x[0]*y[1])]
            for k in (0, 1):
                lo, hi = _bounds(c, bound[k]); assert lo <= got[k] <= hi; count += 1
        Y = c.mpf('1.7')
        for t in (c.mpf(0), c.mpf('.3'), c.mpf(1)):
            ratio = c.exp(4*h*(1-t))
            gr = 4*(1-t)*c.exp(4*c.mpf('.01')*(1-t))
            gm = 4*c.exp(4*c.mpf('.01'))+2/Y*(1+c.mpf('.01')*4*c.exp(4*c.mpf('.01')))
            assert abs(ratio-1)/h <= gr+c.mpf('1e-95'); geometry += 1
            assert abs(ratio*(1-2*h/Y)-1)/h <= gm; geometry += 1
    return dict(passed=True, ordinary_C0_C1_hb_quadratic_error_products=count,
                independent_true_radius_and_length_error_conditions=geometry)


def run():
    start = time.perf_counter()
    source, exponential, products = independent_S1(), independent_exponential_integrals(), independent_error_products()
    print('Original S1, direct macro density, resonant kernels and factored error checks pass', flush=True)
    f = OriginalMacroM15(precision=80)
    audit = mp.mp.clone(); audit.dps = 190
    samples = []; inlet = zeros = errors = units = histories = endpoints = 0
    def check(Z, coordinate):
        nonlocal inlet, zeros, errors, units, histories, endpoints
        p = f.evaluate(Z, coordinate); samples.append(p)
        assert all(p['contract'].values()) and len(p['m15_C0_C1_terms']) == 28
        assert p['m15_C0_C1_terms'][:17] == p['incoming_second_micro_m15_terms']; inlet += 1
        assert p['actual_selected_integer'] == f.inlet(Z)['actual_selected_integer']
        assert p['macro_full_prefix_cover']['cumulative_history_charts'] == ['first_micro', 'second_micro', 'frozen_macro']; histories += 1
        for term in p['m15_C0_C1_terms'][17:]:
            for row in term['coefficient_C0_C1']:
                ends = [audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]
                assert max(abs(x-audit.mpf(row['approximation'])) for x in ends) <= audit.mpf(row['absolute_error_upper']); errors += 1
                if p['exact_empty_increment']:
                    assert ends[0] == ends[1] == audit.mpf(row['approximation']) == 0; zeros += 1
            if 'source_geometry_error' in term['name'] or 'leading_viscous' in term['name']:
                assert term['log_scale_terms'][-1]['source_log_atom_only']; units += 1
        if not p['exact_empty_increment']:
            budget = p['incoming_second_endpoint_error_budget']
            assert budget is not None
            lowlog = audit.make_mpf(tuple(budget['lower_log_radius']['source_enclosure_exact_mpf_tuples'][1]))
            for name, key, power in (('endpoint', 'lower_endpoint', 2), ('viscous', 'lower_viscous', 1)):
                for k in (0, 1):
                    low = budget[key][k]
                    if low['exact_zero']: continue
                    total = p['finite_N_error_logs_in_complete_physical_units'][name][k]
                    assert not total['exact_zero']
                    want = audit.make_mpf(tuple(low['log_absolute_upper']['upper_exact_mpf_tuple']))+power*lowlog
                    got = audit.make_mpf(tuple(total['log_absolute_upper']['upper_exact_mpf_tuple']))
                    # The serialized source arithmetic is interval-directed;
                    # the final atom covers each physical endpoint budget.
                    assert got >= want or abs(got-want) <= abs(got)*audit.mpf('1e-160'); endpoints += 1
            q, _, _ = f.rows.micro._coordinate(coordinate, False)
            lead = f.leading(Z, q)
            assert len(lead['W']) == 3 and _bounds(audit, lead['W'][2])[1] > 0; histories += 1
            assert lead['proxy_R_Qy_constant']
        if coordinate == '1': assert p['macro_geometry']['R100_endpoint_exact']; units += 1
        assert not p['limits']['genuine_temporal_n_ge_2']
        return p
    for Z in ('-.3', '0', '.3'):
        for coordinate in ('0', '1e-20', '1e-18', '.25', '.5', '1'):
            check(Z, coordinate)
        print('Original cumulative m15 reaches physical R=100: Z='+Z, flush=True)
    check('.3', '.875')
    rejected = 0
    for args in [('1', '1'), ('nan', '.5'), ('.3', '-.01'), ('.3', '1.01'), ('.3', 'nan')]:
        try: f.evaluate(*args)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid macro moment query admitted')
    code = dict(f.code_hashes)
    for name in ('scripts/run_original_first_micro_m15.py', 'scripts/run_original_macro_m15.py'):
        code[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family, samples=samples, input_hashes=f.hashes, implementation_hashes=code,
        checks=dict(passed=True, original_source_S1=source, independent_macro_density=exponential,
            original_hb_error_products_and_geometry=products, actual_macro_packets=len(samples),
            exact_complete_second_inlet_replays=inlet, serialized_increment_error_conditions=errors,
            exact_empty_macro_rows=zeros, original_hb_unit_and_R100_conditions=units,
            complete_history_and_resonance_conditions=histories, actual_lower_pressure_endpoint_conditions=endpoints,
            rejected_invalid_queries=rejected), execution_seconds=round(time.perf_counter()-start, 3),
        limits=dict(macro_cumulative_leading_signed_m15_with_complete_factored_errors=True,
            resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
            remaining_active_compact_signed_pressure=False, independent_preRp_signed_seam=False,
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
