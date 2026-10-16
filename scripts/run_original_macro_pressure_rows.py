"""Build original macro pressure rows before cumulative S1 integration."""
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
from openai_ns_reconstruction.original_macro_pressure_rows import OriginalMacroPressureRows
from openai_ns_reconstruction.original_first_micro_pressure_rows import pressure_coefficients

NAME = PREFIX+'compliant_current_original_macro_pressure_rows.json'
MACRO_SOURCE = PREFIX+'compliant_current_original_macro_finite_N.py'


def source_identities():
    """Execute original macro_cell source assignments with independent algebra."""
    tree = ast.parse((DIRECTORY/MACRO_SOURCE).read_text(encoding='utf-8'))
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'macro_cell')
    names = {'modes', 'inv_barphi', 'quotient', 'drives', 'Vy'}
    nodes = [n for n in fn.body if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id in names for t in n.targets)]
    assert len(nodes) == 5
    hb, P2, F2, R, R0, corephi, barq, ell, z = sy.symbols('hb P2 F2 R R0 corephi barq ell Z', positive=True)
    drive = {part: [(j+1)+(j+2)*z+z**(j+2) for j in range(3)]
             for part in ('hydro', 'pressure', 'swirl')}
    def factor(powers):
        assert powers[0] == 1 and powers[3:] == (0, 0)
        return hb*P2**powers[1]*F2**powers[2]
    f = SimpleNamespace(q=barq, multiply=lambda a, b: a*b, jet=lambda x: x,
                        scale=lambda a, b: a*b, add=lambda *rows: sum(rows),
                        factor=factor, drive=drive)
    env = dict(f=f, R0=R0, R=R, phi=corephi*sy.exp(ell),
               core=SimpleNamespace(reciprocal=lambda: 1/corephi),
               view=SimpleNamespace(radial_power=lambda r, j: r**j),
               fields=SimpleNamespace(PARTS=tuple(drive)))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(DIRECTORY/MACRO_SOURCE), 'exec'), env)
    expected = -hb*barq*sy.exp(ell)*sum(
        scale*sum(row*R0**j*R**(power-j) for j, row in enumerate(drive[part]))
        for part, power, scale in (('hydro', 1, 1), ('pressure', 1, P2), ('swirl', 2, F2)))
    for k in (0, 1):
        assert sy.cancel(sy.diff(env['Vy']-expected, z, k)) == 0
    # Pressure recovery consumes the mean ODE before differentiation.
    delta = sy.symbols('delta')
    u, m, uy = [sy.Function(name)(z) for name in ('U', 'M', 'Uy')]
    d, L = 1-z*z, 1-delta*z*z
    q = (2*z*u-(1-delta)*z*m-d*sy.diff(m, z))/L
    qy = (2*z*uy-(1-delta)*z*(u-m)-d*sy.diff(u-m, z))/L
    total = ((1+delta)*z*u+2*z*uy-d*sy.diff(u, z))/L
    for k in (0, 1):
        assert sy.cancel(sy.diff(q+qy-total, z, k)) == 0
    n, by, bphi, byz, bzphi = sy.symbols('N B_y B_phi B_yZ B_Zphi', nonzero=True)
    assert sy.expand((by+n*bphi)/n-(by/n+bphi)) == 0
    assert sy.expand((byz+n*bzphi)/n-(byz/n+bzphi)) == 0
    return dict(passed=True,actual_original_macro_drive_C0_C1_identities=2,
                original_mean_ODE_Q_plus_Qy_C0_C1_identities=2,
                fast_phase_chain_identities=2,
                source_assignments_executed=sorted(names), source=MACRO_SOURCE)


def independent_mean():
    """C0..3 Volterra quadrature with nonzero inlet and evolving source."""
    c = mp.mp.clone(); c.dps = 85
    count = 0
    for text in ('.1', '.7', '4'):
        length = c.mpf(text)
        inlet = lambda z: 2+z+z*z+z**3
        baseline = lambda z: 3-z+z*z/2-z**3/3
        amplitude = lambda z: 1+z+z*z+z**3
        exact = lambda z: (c.exp(-length)*inlet(z)+c.quad(
            lambda y: c.exp(-(length-y))*(baseline(z)+amplitude(z)*(1-c.exp(-y))), [0, length]))
        proxy = lambda z: c.exp(-length)*inlet(z)+(1-c.exp(-length))*baseline(z)
        for k in range(4):
            z = c.mpf('.31')
            err = abs(c.diff(exact, z, k)-c.diff(proxy, z, k))
            bound = (1-c.exp(-length))*abs(c.diff(amplitude, z, k))
            assert 0 < err <= bound
            count += 1
    return count


def run():
    start = time.perf_counter()
    identities, mean = source_identities(), independent_mean()
    print('Original macro drive, pressure identities and independent Z3 Volterra integrals pass', flush=True)
    f = OriginalMacroPressureRows(precision=85)
    audit = mp.mp.clone(); audit.dps = 190
    samples = []; serialized = ordinary = fast = endpoint = history = geometry = high = source = 0
    def check(Z, coordinate):
        nonlocal serialized, ordinary, fast, endpoint, history, geometry, high, source
        p = f.evaluate(Z, coordinate); samples.append(p)
        q, _, _ = f.micro._coordinate(coordinate, False)
        point, truth, geo = f.background(Z, q)
        v = f.inputs(Z, q)
        assert all(p['contract'].values())
        assert p['selected_integer'] == f.micro.integer
        assert p['source_cover']['generic_derivative_chart'] == 'bridge_macro'
        assert p['source_cover']['cumulative_history_charts'] == ['first_micro', 'second_micro', 'frozen_macro']
        assert [len(p['physical_pressure_rows'][k]['signed_background_ordinary_Z'])
                for k in ('U', 'm', 'Q', 'U_y', 'Q_y')] == [3, 4, 3, 2, 2]
        for name, field in p['physical_pressure_rows'].items():
            for k, row in enumerate(field['signed_background_ordinary_Z']):
                actual = truth[name][k]*math.factorial(k)
                assert row['source_enclosure_exact_mpf_tuples'] == [list(t) for t in actual._mpi_]
                ordinary += 1
        groups = [field['signed_background_ordinary_Z'] for field in p['physical_pressure_rows'].values()]
        groups += list(p['normalized_background_ordinary_Z'].values())
        groups += [p['leading_U_y_factored']['coefficient_ordinary_Z']]
        for rows in groups:
            for row in rows:
                ends = [audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]
                assert max(abs(e-audit.mpf(row['approximation'])) for e in ends) <= audit.mpf(row['absolute_error_upper'])
                serialized += 1
        for k, key in enumerate(('B_phi', 'B_Z_phi')):
            row = p['finite_N_axial_radial_error_components']['fast_B_phi'][k]
            want = f.read_log(f.loop['phase_and_first_slow_phase_log_bounds'][key])+f.o.core.lp_iv
            assert row['log_absolute_upper']['lower_exact_mpf_tuple'] == list(want._mpi_[0])
            assert row['log_absolute_upper']['upper_exact_mpf_tuple'] == list(want._mpi_[1])
            fast += 1
        # Full meanZ3 uncertainty survives even if its selected-N error is
        # too large to exponentiate numerically. QZZ must retain it.
        assert len(p['physical_pressure_rows']['m']['finite_N_absolute_error_log_bounds']) == 4
        assert not p['physical_pressure_rows']['m']['finite_N_absolute_error_log_bounds'][3]['exact_zero']
        assert not p['physical_pressure_rows']['Q']['finite_N_absolute_error_log_bounds'][2]['exact_zero']
        high += 2
        ep = p['R_Q_y_endpoint']
        assert ep['signed_background_ordinary_Z'] == p['normalized_background_ordinary_Z']['Qpsi_y']
        for row, value in ((ep['log_physical_factor'], geo['logR']-f.o.core.ll_iv),
                           (ep['finite_N_log_physical_factor'], geo['logR'])):
            assert row['source_enclosure_exact_mpf_tuples'] == [list(t) for t in value._mpi_]
            endpoint += 1
        # Native q=1 is exactly physical R=100, despite enormous rho.
        radius_row = p['geometry']['log_physical_radius']
        lo, hi = [audit.make_mpf(tuple(t)) for t in radius_row['source_enclosure_exact_mpf_tuples']]
        if coordinate == '1':
            assert lo <= audit.log(100) <= hi
            assert p['geometry']['R100_endpoint_exact']; geometry += 2
        if coordinate == '0':
            assert v['truth']['Mpsi'].coefficients == v['incoming']['actual']['Mpsi'].coefficients
            history += 1
        assert not p['limits']['macro_signed_m15_integrated']
        assert all(v['point'][name].order >= 3 and v['truth'][name].order >= 3 for name in ('Psi', 'Mpsi'))
        # Actual radial source is exported before hb, rather than replacing
        # the unknown original width by its numerical error cap.
        assert v['source_Uy']['truth'].order == 5
        assert set(v['source_Uy']['parts']) == {'hydro', 'pressure', 'swirl'}
        assert any(_bounds(audit, abs(t))[1] > 0 for t in v['source_Uy']['truth'].coefficients)
        source += 3
        return p
    for Z in ('-.3', '0', '.3'):
        for coordinate in ('0', '1e-20', '1e-18', '.25', '.5', '1'):
            check(Z, coordinate)
        a = f.inputs(Z, f.iv.mpf(0))['point']['Mpsi'][0]
        b = f.inputs(Z, f.iv.mpf('1e-18'))['point']['Mpsi'][0]
        assert a != b; history += 1
        print('Original macro QZZ, Uy/Qy and actual R100 endpoint: Z='+Z, flush=True)
    check('.3', '.875')
    # Isolate the genuine third-mean channel with the exact pressure
    # coefficient. This cannot pass with a padded zero M_ZZZ row.
    errors = {name: [None]*4 for name in ('U', 'm', 'U_y')}
    errors['m'][3] = f.iv.mpf(0)
    errors['U_Z'], errors['m_Z'] = errors['U'][1:], errors['m'][1:]
    cc = pressure_coefficients(f.J, f.iv, '.3', f.o.core.delta_iv, 2)
    isolated = f.combine(cc['Q'], errors, 2)
    assert isolated[0] is None and isolated[1] is None and isolated[2] is not None; high += 1
    rejected = 0
    for Z, q in (('1', '1'), ('nan', '.5'), ('.3', '-.01'), ('.3', '1.01'), ('.3', 'nan')):
        try: f.evaluate(Z, q)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid macro pressure query admitted')
    try: f.evaluate('.3', '1', selected_sc_multiple=True)
    except ValueError: rejected += 1
    else: raise AssertionError('Selected-sc coordinate admitted as macro')
    for name in (MACRO_SOURCE,):
        digest = hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()
        assert f.hashes.get(name, digest) == digest
        f.hashes[name] = digest
    code = dict(f.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family, samples=samples, input_hashes=f.hashes, implementation_hashes=code,
        checks=dict(passed=True, original_source_identities=identities, independent_mean_Z0_Z3_integrals=mean,
                    actual_macro_pressure_packets=len(samples), ordinary_physical_rows=ordinary,
                    serialized_error_conditions=serialized, original_fast_phase_conditions=fast,
                    finite_R_Qy_unit_conditions=endpoint, nonzero_inlet_and_mean_relaxation_conditions=history,
                    exact_R100_log_geometry_conditions=geometry, genuine_mean_Z3_and_QZZ_conditions=high,
                    actual_macro_modal_drive_conditions=source, rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start, 3),
        limits=dict(actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors=True,
                    resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
                    macro_signed_m15_integrated=False, independent_preRp_signed_seam=False,
                    global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
                    stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))


def inspect(path):
    d = json.loads(path.read_text(encoding='utf-8')); assert d['checks']['passed']
    for n, h in d['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest() == h, n
    for n, h in d['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest() == h, n
    return d


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write-report', action='store_true'); p.add_argument('--inspect-report', action='store_true')
    args = p.parse_args(); path = DIRECTORY/NAME
    d = inspect(path) if args.inspect_report else run()
    if args.write_report:
        path.write_text(json.dumps(d, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(checks=d['checks'], execution_seconds=d['execution_seconds']), indent=2), flush=True)
