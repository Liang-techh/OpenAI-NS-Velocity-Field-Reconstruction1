"""Build actual first-switch pressure rows with native hbS controls."""
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
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_first_switch_pressure_rows import OriginalFirstSwitchPressureRows

NAME = PREFIX+'compliant_current_original_first_switch_pressure_rows.json'
SOURCE = PREFIX+'compliant_current_original_prescribed_R100_Rm_C6_functions.py'


def source_controls():
    """Execute original source bodies using independent symbolic algebra."""
    tree = ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('full_drive', 'first_switch_fields')]
    assert len(nodes) == 2
    s, z, h, F0, Pstar = sy.symbols('s Z hbS F0 Pstar', positive=True)
    radius = 100*sy.exp(h*s)
    phi0, barphi = 2+z+z*z, 3+z*z
    known = dict(D_over_R=1+z/3, drive_hydro=1+z*z,
                 drive_pressure=2-z/2, drive_swirl=1+z/5)
    g = SimpleNamespace(symbol=sy.Symbol, constant=sy.Rational, one=sy.Integer(1), zero=sy.Integer(0),
                        mul=lambda *v: sy.prod(v), sub=lambda a, b: a-b, neg=lambda v: -v)
    a = SimpleNamespace(g=g, scale=lambda a, b: a*sy.sympify(b), mul=lambda *v: sy.prod(v),
        add=lambda *v: sum(v), exp=sy.exp, trim=lambda v: v, div=lambda a, b, reason: a/b,
        at=lambda value, variable, argument: value.subs(variable, argument),
        integral=lambda value, variable, lower, upper, **kwargs: sy.Integral(value, (variable, lower, upper)))
    sigma = lambda graph, phase: sy.Function('sigma')(phase)
    env = dict(sigma=sigma)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(DIRECTORY/SOURCE), 'exec'), env)
    result = env['first_switch_fields'](a, known, dict(phi=phi0, U=4*z+1), s, h, radius, F0, Pstar, barphi)
    drive = radius*known['drive_hydro']+radius*Pstar**2*known['drive_pressure']+radius**2*F0**2*known['drive_swirl']
    expressions = [result['U_y']+h*(1-sigma(g, s))*result['phi']/barphi*drive,
                   result['phi_y']+h*radius*known['D_over_R']*result['phi']/2,
                   result['angular_phase_density']+h*h*radius*known['D_over_R']/2]
    count = 0
    for expr in expressions:
        for k in (0, 1):
            assert sy.cancel(sy.diff(expr, z, k)) == 0; count += 1
    assert sy.simplify(sy.diff(result['U'], s)-h*result['U_y']) == 0
    return dict(passed=True, original_first_switch_source_C0_C1_identities=count,
                actual_phase_velocity_ODE_with_two_native_width_factors=True,
                native_hbS_is_parameter_not_power=True, source_bodies_executed=['full_drive', 'first_switch_fields'], source=SOURCE)


def independent_mean():
    c = mp.mp.clone(); c.dps = 85
    count = 0
    for h in (c.mpf('.003'), c.mpf('.02')):
        for s in (c.mpf('.2'), c.mpf('.8')):
            initial = lambda z: 2+z+z*z+z**3
            psi0 = lambda z: 3-z+z*z/2-z**3/3
            amplitude = lambda z: 1+z+z*z+z**3
            exact = lambda z: c.exp(-h*s)*initial(z)+h*c.quad(
                lambda t: c.exp(-h*(s-t))*(psi0(z)+h*h*amplitude(z)*t), [0, s])
            for k in range(4):
                z = c.mpf('.31')
                error = abs(c.diff(exact, z, k)-c.diff(initial, z, k))
                bound = h*s*(abs(c.diff(lambda z: psi0(z)-initial(z), z, k))
                             +h*h*s*abs(c.diff(amplitude, z, k)))
                assert 0 < error <= bound; count += 1
    return count


def run():
    start = time.perf_counter(); controls, mean = source_controls(), independent_mean()
    print('Actual first-switch controls, two-width source identity and independent meanZ0..3 integrals pass', flush=True)
    f = OriginalFirstSwitchPressureRows(precision=85)
    audit = mp.mp.clone(); audit.dps = 190
    samples = []; ordinary = serialized = fast = inlet = zeros = width = endpoint = 0
    def check(Z, coordinate):
        nonlocal ordinary, serialized, fast, inlet, zeros, width, endpoint
        p = f.evaluate(Z, coordinate); samples.append(p)
        s, _, _ = f.micro._coordinate(coordinate, False)
        point, truth, geo = f.background(Z, s); v = f.inputs(Z, s)
        assert all(p['contract'].values())
        assert f.evaluate(Z, coordinate) is p
        assert geo['length_is_error_cover_only'] and geo['theta_mass_are_error_covers_only']
        assert p['selected_integer'] == f.micro.integer
        assert p['source_cover']['generic_derivative_chart'] == 'switch_first'
        assert p['source_cover']['cumulative_history_charts'] == ['first_micro', 'second_micro', 'frozen_macro', 'first_switch']
        assert f.macro_rows.ctx is f.ctx and f.macro_rows.iv is f.iv
        for name, field in p['physical_pressure_rows'].items():
            for k, row in enumerate(field['signed_background_ordinary_Z']):
                value = truth[name][k]*math.factorial(k)
                assert row['source_enclosure_exact_mpf_tuples'] == [list(t) for t in value._mpi_]; ordinary += 1
        groups = [field['signed_background_ordinary_Z'] for field in p['physical_pressure_rows'].values()]
        groups += list(p['normalized_background_ordinary_Z'].values())
        groups += [p['leading_U_y_factored']['coefficient_ordinary_Z']]
        for rows in groups:
            for row in rows:
                ends = [audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]
                assert max(abs(e-audit.mpf(row['approximation'])) for e in ends) <= audit.mpf(row['absolute_error_upper']); serialized += 1
        evidence = p['first_switch_source_control_evidence']
        assert evidence['angular_weight'] == '1' and evidence['axial_weight'] == '1-sigma(s)'
        assert evidence['original_hb_power_angular'] == evidence['original_hb_power_axial'] == 2; width += 1
        atom = p['geometry']['original_log_width']
        assert atom['original_native_parameter_name'] == 'hbS' and atom['bridge_hbB_graph_node_not_aliased']; width += 1
        for gate in ('shared_upstream_is_same_object', 'bridge_logh_is_actual_upstream_object',
                     'switch_width_is_exact_source_reboxing', 'bridge_hbB_switch_hbS_kept_distinct'):
            assert p['geometry']['width_source_identity'][gate]; width += 1
        assert not p['physical_pressure_rows']['m']['finite_N_absolute_error_log_bounds'][3]['exact_zero']
        assert not p['physical_pressure_rows']['Q']['finite_N_absolute_error_log_bounds'][2]['exact_zero']
        for k, key in enumerate(('B_phi', 'B_Z_phi')):
            got = p['finite_N_axial_radial_error_components']['fast_B_phi'][k]
            want = f.read_log(f.loop['phase_and_first_slow_phase_log_bounds'][key])+f.o.core.lp_iv
            assert got['log_absolute_upper']['lower_exact_mpf_tuple'] == list(want._mpi_[0])
            assert got['log_absolute_upper']['upper_exact_mpf_tuple'] == list(want._mpi_[1]); fast += 1
        ep = p['R_Q_y_endpoint']
        assert ep['signed_background_ordinary_Z'] == p['normalized_background_ordinary_Z']['Qpsi_y']
        for row, actual in ((ep['log_physical_factor'], geo['logR']-f.o.core.ll_iv),
                            (ep['finite_N_log_physical_factor'], geo['logR'])):
            assert row['source_enclosure_exact_mpf_tuples'] == [list(t) for t in actual._mpi_]; endpoint += 1
        if coordinate == '0':
            upstream = f.macro_rows.evaluate(Z, '1')
            for name in ('Psi', 'Mpsi', 'Qpsi'):
                assert p['normalized_background_ordinary_Z'][name] == upstream['normalized_background_ordinary_Z'][name]; inlet += 1
            assert p['leading_U_y_factored']['coefficient_ordinary_Z'] == upstream['leading_U_y_factored']['coefficient_ordinary_Z']; inlet += 1
            assert p['geometry']['R100_inlet_exact']; inlet += 1
        if coordinate == '1':
            assert p['geometry']['R1_endpoint_exact']
            assert v['source_Uy']['first_switch_complement']._mpi_ == f.iv.mpf(0)._mpi_
            for k in (0, 1):
                row = p['leading_U_y_factored']['coefficient_ordinary_Z'][k]
                assert audit.mpf(row['approximation']) == 0
                assert all(audit.make_mpf(tuple(t)) == 0 for t in row['source_enclosure_exact_mpf_tuples']); zeros += 1
            assert not p['physical_pressure_rows']['U_y']['finite_N_absolute_error_log_bounds'][0]['exact_zero']
        assert not p['limits']['first_switch_signed_m15_integrated']
        return p
    for Z in ('-.3', '0', '.3'):
        for coordinate in ('0', '.125', '.25', '.5', '.75', '1'):
            check(Z, coordinate)
        print('Original first-switch pressure/control/endpoint rows through R1: Z='+Z, flush=True)
    check('.3', '.875')
    rejected = 0
    for args in [('1', '1'), ('nan', '.5'), ('.3', '-.01'), ('.3', '1.01'), ('.3', 'nan')]:
        try: f.evaluate(*args)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid original first-switch request admitted')
    try: f.evaluate('.3', '.5', selected_sc_multiple=True)
    except ValueError: rejected += 1
    else: raise AssertionError('Selected-sc coordinate admitted as first switch')
    digest = hashlib.sha256((DIRECTORY/SOURCE).read_bytes()).hexdigest()
    assert f.hashes.get(SOURCE, digest) == digest; f.hashes[SOURCE] = digest
    code = dict(f.code_hashes); code[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family, samples=samples, input_hashes=f.hashes, implementation_hashes=code,
        checks=dict(passed=True, original_source_controls=controls, independent_mean_Z0_Z3_integrals=mean,
            actual_first_switch_pressure_packets=len(samples), ordinary_physical_rows=ordinary,
            serialized_error_conditions=serialized, original_fast_phase_conditions=fast,
            exact_macro_inlet_and_factored_derivative_conditions=inlet, exact_terminal_unmodified_Uy_rows=zeros,
            original_hbS_width_identity_and_control_conditions=width, finite_R_Qy_unit_conditions=endpoint,
            rejected_invalid_queries=rejected), execution_seconds=round(time.perf_counter()-start, 3),
        limits=dict(actual_first_switch_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors=True,
            resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
            first_switch_signed_m15_integrated=False, independent_preRp_signed_seam=False,
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
