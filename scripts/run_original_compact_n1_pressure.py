"""Original compact n1 pressure, incoming histories and physical order1."""
import argparse
import ast
import hashlib
import inspect as inspectlib
import json
from pathlib import Path
import sys
import textwrap
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'experiments/root_st073')]
DIRECTORY = ROOT / 'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_compact_n1_pressure import (
    OriginalReservedN1PressureNumeric, OriginalCompactN1PressureVelocity, weighted_beta_prefix)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
from run_original_second_micro_numeric import overlap

NAME = PREFIX + 'compliant_current_original_compact_n1_pressure_numeric.json'


def independent_integrals():
    c = mp.mp.clone(); c.dps = 75
    iv = MPIntervalContext(); iv.dps = 90
    engine = SimpleNamespace(ctx=c, iv=iv, J=J, ell=iv.log(2)/400,
                             cells=256, beta_jets=beta_jets, _pressure_integrals={})
    count = 0; maximum = c.mpf(0)
    for mu in ('0', '.001', '.4'):
        rate = -iv.mpf('.5') - iv.mpf(mu)
        for endpoint in ('-1', '-.75', '0', '.75', '1'):
            got = weighted_beta_prefix(engine, c.mpf(endpoint), rate)
            truth = c.quad(lambda t: c.exp(-1/(1-t*t) -
                           (c.mpf('.5')+c.mpf(mu))*c.log(2)*t/400) if abs(t)<1 else c.mpf(0),
                           [-1, c.mpf(endpoint)])
            lo, hi = _bounds(c, got)
            assert lo-c.mpf('1e-65') <= truth <= hi+c.mpf('1e-65'), (mu, endpoint)
            maximum = max(maximum, hi-lo); count += 1
    return count, c.nstr(maximum, 16)


def independent_source_algebra():
    import sympy as sy
    import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as original
    Rc, lp, lu = sy.symbols('logRc logP logUp', real=True)
    lm = -4*lp-sy.log(1000)
    ctx = SimpleNamespace(mpf=lambda v: sy.Rational(str(v)), log=sy.log, ln=sy.log,
                          exp=sy.exp, sqrt=sy.sqrt, pi=sy.pi)
    namespace = dict(original.pressure_unit_log.__globals__)
    namespace['interval'] = lambda c, v: v
    tree = ast.parse(textwrap.dedent(inspectlib.getsource(original.pressure_unit_log)))
    exec(compile(tree, '<original m13 pressure unit>', 'exec'), namespace)
    target = namespace['pressure_unit_log'](ctx, dict(logP=lp, logUp=lu, logRc=Rc), 1)
    f = object.__new__(OriginalReservedN1PressureNumeric)
    f.iv = ctx; f.reserved = SimpleNamespace(logUp=lu)
    f.background = SimpleNamespace(patch=SimpleNamespace(leading=SimpleNamespace(defect=SimpleNamespace(logP=lp))),
                                   outer=SimpleNamespace(logmui=lm))
    f._row = lambda v: v; f._rows = lambda v: v
    actual = sum(f.pressure_unit_atoms(dict(coordinate=dict(global_y_log_terms=[Rc-sy.log(110), 0, 0, 0, sy.Symbol('s')]))))
    assert sy.simplify(sy.expand(actual-target)) == 0
    R, R0, k = sy.symbols('R Rc k', positive=True)
    Q = -k*R0/R
    # The unweighted reduced identity has only B(R Q^2)/2; its
    # viscous endpoint Q+R Q_R vanishes, including both finite endpoints.
    prefix = -(R*Q**2 - R0*k**2)/4
    expected = R0*k**2*(1-R0/R)/4
    assert sy.simplify(prefix-expected) == 0
    assert sy.simplify(Q+R*sy.diff(Q, R)) == 0
    assert sy.simplify(sy.diff(prefix, R)-Q**2/4) == 0
    return dict(fixed_m13_original_unit_identity=True, exact_partial_baseline_identity=True,
                both_finite_viscous_endpoints_exact_zero_for_baseline=True)


def run():
    start = time.perf_counter()
    count, width = independent_integrals(); algebra = independent_source_algebra()
    print('Independent original pressure-unit/baseline algebra and weighted beta integrals pass', flush=True)
    f = OriginalCompactN1PressureVelocity(); c, iv = f.ctx, f.iv
    samples = []; ownership = finite = products = physical = support = seams = 0
    radii = [('inlet', '1'), ('before_swirl', c.nstr(c.exp(c.log(2)*c.mpf('0.31')), 95))]
    for j in (13, 14, 15):
        for shift, label in ((c.mpf(0), 'center'), (c.mpf(1)/800, 'partial')):
            radii.append((label+str(j), c.nstr(c.exp(c.log(2)*(c.mpf(j)/40+shift)), 95)))
    radii.extend([('suffix', '1.9'), ('terminal', '2')])
    for Z in ('-.3', '0', '.3'):
        for label, x in radii:
            p = f.evaluate(Z, x, log_tau='-2', azimuth='.25')
            p['sample_label'] = label; samples.append(p)
            n = p['reserved_n1_pressure']; terms = {v['name']: v for v in n['P1_C0_C1_terms']}
            assert set(terms) == {'inner', 'lower', 'active', 'compact_finite_N', 'compact_exact_leading', 'own_swirl'}
            assert all(n['contract'].values()); ownership += 1
            assert p['selected_integer'] == f.numeric._inlets[Z]['selected_integer']; ownership += 1
            assert n['cumulative_true_function'] == f.pressure.binding['cumulative']['m13']; ownership += 1
            assert p['physical']['leading_pressure_includes_same_P0_once']; ownership += 1
            for order in n['actual_partial_finite_N_normalized_log_bounds']:
                for row in order:
                    assert _bounds(c, f._restore(row))[1] < -1000; finite += 1
            baseline = terms['compact_exact_leading']['coefficient_C0_C1']
            lo, hi = _bounds(c, f._restore(baseline[0]))
            truth = 1-1/c.mpf(x)
            assert lo-c.mpf('1e-75') <= truth <= hi+c.mpf('1e-75'); products += 1
            assert _bounds(c, f._restore(baseline[1])) == (0, 0); products += 1
            for row in n['swirl_partial_integral_factors']:
                assert _bounds(c, f._restore(row))[0] >= 0; support += 1
            if label in ('inlet', 'before_swirl'):
                for row in terms['own_swirl']['coefficient_C0_C1']:
                    assert _bounds(c, f._restore(row)) == (0, 0); support += 1
            if label == 'inlet':
                for name in ('compact_finite_N', 'compact_exact_leading', 'own_swirl'):
                    for row in terms[name]['coefficient_C0_C1']:
                        assert _bounds(c, f._restore(row)) == (0, 0); seams += 1
                for name in ('inner', 'lower', 'active'):
                    for k in (0, 1):
                        got = f._restore(terms[name]['coefficient_C0_C1'][k])
                        assert _bounds(c, got)[0] < _bounds(c, got)[1]; seams += 1
            assert p['physical']['full_n1_pressure_included']; physical += 1
            assert p['physical']['pressure_orders_included'] == [0, 1]; physical += 1
            assert p['physical']['n1_pressure_time_weight'] == 'lambda^(-2)'; physical += 1
            assert not p['physical']['n1_pressure_extra_Lambda_squared_factor_applied']; physical += 1
            assert not p['physical']['useful_combined_pressure_relative_error_certified']; physical += 1
            assert not n['limits']['terminal_pressure_zero_at_2Rc']; physical += 1
            assert len(p['physical']['pressure']) == 7; physical += 1
        print('Actual compact P1 assembled with all incoming histories: Z='+Z, flush=True)
    aliases = 0
    for x in ('1.0', '1.000e0'):
        p = f.pressure.evaluate('.3', x)
        for term in p['P1_C0_C1_terms']:
            if term['name'] in ('compact_finite_N', 'compact_exact_leading', 'own_swirl'):
                for row in term['coefficient_C0_C1']:
                    assert _bounds(c, f._restore(row)) == (0, 0); aliases += 1
    scaled = 0
    for Z in ('-.3', '0', '.3'):
        _, x = next(v for v in radii if v[0] == 'center14')
        p = f.evaluate(Z, x, log_tau=None, scaled_time='1e18', azimuth='1.1')
        p['sample_label'] = 'scaled_center14'; samples.append(p)
        normal = next(v for v in samples if v['Z'] == Z and v['sample_label'] == 'center14')
        for a, b in zip(p['physical']['pressure'], normal['physical']['pressure']):
            overlap(a['coefficient'], b['coefficient']); scaled += 1
        assert p['physical']['time']['delta_cancelled_before_time_comparison']; scaled += 1
    code = dict(f.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family, samples=samples, input_hashes=f.hashes, implementation_hashes=code,
                checks=dict(passed=True, independent_weighted_beta_integrals=count, max_tested_integral_width=width,
                            independent_source_algebra=algebra, actual_compact_P1_physical_packets=len(samples),
                            source_and_incoming_history_conditions=ownership, actual_partial_finite_N_bounds=finite,
                            exact_nonzero_baseline_conditions=products, support_and_positive_integral_conditions=support,
                            inlet_pressure_history_conditions=seams, physical_pressure_order_and_unit_conditions=physical,
                            scaled_time_conditions=scaled, exact_inlet_alias_conditions=aliases), execution_seconds=round(time.perf_counter()-start, 3),
                limits=dict(compact_velocity_and_pressure_temporal_orders=[0, 1], P1_ordinary_Z1=True,
                            resolved_signed_function_centers=False, relative_pressure_cancellation_accuracy=False,
                            independent_signed_terminal_pressure_residual=False, pressure_radial_derivatives=False,
                            Q1_Z=False, beyond_2Rc=False, global_n1_joins=False, common_global_N=False,
                            exact_heat_exterior=False, stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))


def inspect(path):
    data = json.loads(path.read_text(encoding='utf-8'))
    assert data['checks']['passed']
    for name, digest in data['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest, name
    for name, digest in data['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true')
    parser.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args(); path = DIRECTORY/NAME
    data = inspect(path) if args.inspect_report else run()
    if args.write_report:
        path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(checks=data['checks'], execution_seconds=data['execution_seconds']), indent=2), flush=True)
