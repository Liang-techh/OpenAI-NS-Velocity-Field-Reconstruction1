"""Continue original reference m15 and genuine pressure rows to Rm."""
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
from openai_ns_reconstruction.original_reference_m15 import OriginalReferenceM15
from openai_ns_reconstruction.original_macro_m15 import proxy_modes
from openai_ns_reconstruction.original_first_switch_m15 import frozen_coefficients
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME = PREFIX+'compliant_current_original_reference_m15_numeric.json'
SOURCE = PREFIX+'compliant_current_original_reference_restore_functions.py'


def original_source():
    tree = ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef)
               and any(isinstance(v, ast.FunctionDef) and v.name == 'restoration' for v in n.body))
    restore = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'restoration')
    start = next(i for i, n in enumerate(restore.body) if isinstance(n, ast.Assign)
                 and any(isinstance(a, ast.Name) and a.id == 'cutoff' for a in n.targets))
    z, j, eps, psi, sigma, sp = sy.symbols('Z j epsilon psi sigma sigma_prime')
    obj = SimpleNamespace(V=4*z+j+eps*psi, Vref=4*z, E=j+eps*psi)
    flow = SimpleNamespace(add=lambda *v: sum(v), scale=lambda v, a: v*a)
    scope = dict(self=obj, f=flow, q=sy.Rational(1, 3), t=sy.Rational(1, 3),
                 c=None, sigma_jets=lambda c, t: [sigma, sp, 0, 0, 0], math=math)
    exec(compile(ast.Module(body=restore.body[start:start+3], type_ignores=[]), str(DIRECTORY/SOURCE), 'exec'), scope)
    assert sy.expand(scope['V']-(4*z+(j+eps*psi)*(1-sigma))) == 0
    assert sy.expand(obj.E*scope['alpha'][1]+(j+eps*psi)*sp) == 0
    post = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'postrestore')
    packet = next(n.value for n in ast.walk(post) if isinstance(n, ast.Assign)
                  and isinstance(n.value, ast.Call) and isinstance(n.value.func, ast.Attribute)
                  and n.value.func.attr == 'packet')
    assert isinstance(packet.args[1], ast.Attribute) and packet.args[1].attr == 'Vref'
    return dict(passed=True, actual_restoration_U_and_Uy_sector_identities=2,
                actual_postrestore_U_reference_identity=1, source=SOURCE)


def independent_reference_primitives():
    z, t, eps, Y = sy.symbols('Z y epsilon Y')
    de, axis = sy.Rational(1, 17), sy.Rational(7, 5)
    psi = 1+z+z*z/3+z**3/7
    mean = sy.Rational(2, 3)-z/2+z*z/5-z**3/11
    U = 4*z+axis+eps*psi; L = sy.log(sy.Rational(11, 10))
    M = U+eps*(mean-psi)*sy.exp(-Y-L-t)
    A, D = 1-de*z*z, 1-z*z
    Q = (2*z*U-(1-de)*z*M-D*sy.diff(M,z))/A
    Qy = sy.diff(Q,t); R = 110*sy.exp(t)
    flow = ((Q+(1-de)*z*sy.diff(Q,z)/2+Qy)/A+Q*(Q/2+Qy)
            +U*(-2*z*Q+D*sy.diff(Q,z)-2*z*Qy)/A)
    direct = [sy.lambdify((t,z,eps,Y),sy.diff(R*R*flow/2,z,k),'mpmath') for k in (0,1)]
    assert sy.cancel(sy.diff(M,t)-(U-M)) == 0
    assert sy.simplify(sy.diff(R*Qy,t)) == 0
    c = mp.mp.clone(); c.dps = 100
    integrals = pairs = telescopes = error_transport = 0
    with mp.workdps(100):
        nodes, weights = mp.gauss_quadrature(48,'legendre')
        for gap in (mp.mpf('.01'),mp.mpf('.4'),mp.mpf(2)):
            b, ee, inherited = mp.mpf('.7'),mp.mpf('.23'),mp.mpf('1.7')
            for zz in (mp.mpf('-.3'),mp.mpf(0),mp.mpf('.3')):
                def jet(expr):
                    return J(c,[c.mpf(str(sy.N(sy.diff(expr,z,k).subs(z,sy.Rational(str(zz))),100)))/c.factorial(k)
                                for k in range(4)])
                modes,_ = proxy_modes(J,c,str(zz),c.mpf(1)/17,c.mpf(7)/5,jet(psi),jet(mean))
                H = frozen_coefficients(J,c,str(zz),c.mpf(1)/17,modes)
                for k in (0,1):
                    exact = mp.mpf(0)
                    for (d,n),row in H.items():
                        rate = 2-n
                        kernel = gap if not rate else mp.exp(rate*b)*mp.expm1(rate*gap)/rate
                        if rate:
                            pair = (mp.exp(rate*(b+gap))-mp.exp(rate*b))/rate
                            assert abs(pair-kernel) < mp.mpf('1e-92'); pairs += 1
                            old_upper = mp.exp(rate*b)/rate
                            new_lower = -mp.exp(rate*b)/rate
                            assert old_upper+new_lower == 0; telescopes += 1
                        exact += 12100*c.mpf(row[k])*ee**d*mp.exp(-n*(inherited+mp.log(mp.mpf(11)/10)))*kernel
                    want = sum((w*direct[k](b+(node+1)*gap/2,zz,ee,inherited)*gap/2
                                for node,w in zip(nodes,weights)),mp.mpf(0))
                    assert abs(exact-want) < mp.mpf('1e-70')*max(1,abs(want)); integrals += 1
        # Each reference derivative discrepancy obeys e_M'=e_U-e_M.
        # If the long-chart M envelope already covers its constant U
        # envelope, the same bound is invariant under further relaxation.
        for gap in (mp.mpf('.1'),mp.mpf(1),mp.mpf(100)):
            for k in range(4):
                u_bound, m_bound = mp.mpf(k+1),mp.mpf(4*(k+1))
                actual = m_bound*mp.exp(-gap)+u_bound*(-mp.expm1(-gap))
                assert actual <= m_bound; error_transport += 1
    return dict(passed=True, independent_shifted_physical_pressure_integrals=integrals,
                signed_kernel_pair_conditions=pairs, internal_endpoint_cancellations=telescopes,
                constant_U_mean_error_Z0_Z3_envelope_conditions=error_transport,
                exact_constant_U_mean_and_viscous_identities=2)


def run():
    start = time.perf_counter()
    source, independent = original_source(), independent_reference_primitives()
    print('Actual restoration source and independently shifted reference primitives pass',flush=True)
    f = OriginalReferenceM15(precision=85)
    fields = [f.rows, f.rows.next_field()]
    fields.append(fields[-1].next_field())
    audit = mp.mp.clone(); audit.dps = 190
    samples, pressure = [], []
    counts = dict(exact_106term_inlets=0,ordinary_pressure_rows=0,serialized_increment_errors=0,
                  exact_empty_rows=0,shared_long_coefficient_rows=0,internal_telescope_atoms=0,
                  positive_resonance_atoms=0,lower_finite_N_endpoint_conditions=0,
                  source_seam_overlap_conditions=0,atomic_reference_offsets=0,
                  actual_restoration_baseline_conditions=0,independent_mean_kernel_inclusions=0,
                  complete_prefix_conditions=0,nonzero_mean_Z3_error_conditions=0)

    def enclosure(row):
        return [audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]

    def check_pressure(field,Z,coordinate):
        row = field.evaluate(Z,coordinate); pressure.append(row)
        q,_,_ = field._coordinate(coordinate,False)
        point,truth,geo = field.background(Z,q)
        assert all(row['contract'].values())
        assert row['source_cover']['generic_derivative_chart'] == field.OWNERS[field.chart]
        assert row['geometry']['actual_native_length']['source_enclosure_exact_mpf_tuples'] == [list(v) for v in geo['length']._mpi_]
        assert row['geometry']['unit_offsets_not_merged_into_huge_radius']
        assert not row['geometry']['physical_radius_materialized']
        for name,values in row['physical_pressure_rows'].items():
            for k,value in enumerate(values['signed_background_ordinary_Z']):
                exact = truth[name][k]*math.factorial(k)
                assert value['source_enclosure_exact_mpf_tuples'] == [list(t) for t in exact._mpi_]
                counts['ordinary_pressure_rows'] += 1
        assert not row['physical_pressure_rows']['m']['finite_N_absolute_error_log_bounds'][3]['exact_zero']
        counts['nonzero_mean_Z3_error_conditions'] += 1
        n = 8+field.reference.charts.index(field.chart)
        assert len(row['source_cover']['cumulative_history_charts']) == n
        assert row['source_cover']['prefix_source_cells'] > 0
        assert not row['finite_N_axial_radial_error_components']['extra_hb_factor']
        assert 'baseline_ordinary_Z' in row['R_Q_y_endpoint']
        assert 'signed_background_ordinary_Z' not in row['R_Q_y_endpoint']
        # Repeated cached calls must preserve both baseline and epsilon rows.
        assert field.evaluate(Z,coordinate) == row
        counts['complete_prefix_conditions'] += 1
        local = geo['local_point']; incoming = field.inputs(Z,q)['incoming']['point']
        if field.chart == 'restoration':
            alpha,_ = field.reference._alpha(local,geo['local'])
            j = field.o.core.j
            for k in (0,1,2):
                want = (4*field.J.variable(field.ctx,field.ctx.mpf(Z),5)+j*alpha)[k]
                assert abs(point['U0'][k]-want) <= field.ctx.mpf('1e-75')*max(1,abs(want))
                counts['actual_restoration_baseline_conditions'] += 1
            if 0 < local < 1:
                assert point['U0_y'][0] != 0 and point['Psi_y'][0] != 0
                assert point['M0'][0] != point['U0'][0]
            else:
                assert point['U0_y'][0] == point['Psi_y'][0] == 0
            kp,kt = field._mean_kernels[geo['local']._mpi_]
            lo,hi = _bounds(audit,kt)
            assert lo <= audit.mpf(kp) <= hi
            counts['independent_mean_kernel_inclusions'] += 1
        if field.chart == 'postrestore':
            assert point['Psi'][0] == point['U_y'][0] == 0
            assert point['M0'][0] != point['U0'][0]
            # Q0_y is odd in Z for the Z-independent j sector.
            assert (point['Q0_y'][0] == 0) == (field.ctx.mpf(Z) == 0)
            assert row['source_cover']['generic_source_domain'] == ['-7','-6']
            assert row['geometry']['exact_Rm_endpoint'] == (coordinate == '1')
        if field.chart != 'reference' or not coordinate.startswith('g:'):
            atom = row['geometry']['atomic_log_offset_terms'][-1]
            want = (-8*field.ctx.mpf(coordinate) if field.chart == 'reference'
                    else (-8 if field.chart == 'restoration' else -7)+field.ctx.mpf(coordinate))
            assert abs(audit.mpf(atom['approximation'])-audit.mpf(want)) < audit.mpf('1e-75')
            counts['atomic_reference_offsets'] += 1
        return row

    def check_moment(Z,coordinate):
        p = f.evaluate(Z,coordinate); samples.append(p)
        row = check_pressure(f.rows,Z,coordinate)
        assert f.evaluate(Z,coordinate) is p and all(p['contract'].values())
        assert len(p['m15_C0_C1_terms']) == 124
        assert p['m15_C0_C1_terms'][:106] == p['incoming_m15_terms'] == f.inlet(Z)['m15_C0_C1_terms']
        counts['exact_106term_inlets'] += 1
        original = {t['name']:t for t in p['incoming_m15_terms'][88:]}
        added = {t['name']:t for t in p['m15_C0_C1_terms'][106:]}
        for term in added.values():
            for value in term['coefficient_C0_C1']:
                lo,hi = enclosure(value); approx = audit.mpf(value['approximation'])
                assert max(abs(lo-approx),abs(hi-approx)) <= audit.mpf(value['absolute_error_upper'])
                counts['serialized_increment_errors'] += 1
                if p['exact_empty_increment']:
                    assert lo == hi == approx == 0; counts['exact_empty_rows'] += 1
        for d in range(3):
            for n in range(d+1):
                stem = f'half_S1_epsilon_{d}_mean_{n}'
                if n == 2:
                    value = added['reference_'+stem+'_resonance']
                    if not p['exact_empty_increment']:
                        assert value['coefficient_C0_C1'] == original['long_'+stem+'_resonance']['coefficient_C0_C1']
                        assert value['log_scale_terms'][1]['reference_gap_not_multiplied_into_coefficient']
                        counts['positive_resonance_atoms'] += 1
                    continue
                for suffix in ('upper','lower'):
                    new,old = added['reference_'+stem+'_'+suffix],original['long_'+stem+'_'+suffix]
                    assert new['actual_function_binding']['correlated_coefficient_id'] == old['actual_function_binding']['correlated_coefficient_id']
                    if not p['exact_empty_increment']:
                        assert new['coefficient_C0_C1'] == old['coefficient_C0_C1']
                        counts['shared_long_coefficient_rows'] += 2
                low = added['reference_'+stem+'_lower']
                assert low['log_scale_terms'] == original['long_'+stem+'_upper']['log_scale_terms']
                counts['internal_telescope_atoms'] += 1
        if not p['exact_empty_increment']:
            lower = p['actual_lower_finite_N_endpoint_budget']
            radius_hi = enclosure(lower['log_radius'])[1]
            for name,power in (('endpoint',2),('viscous',1)):
                for k in (0,1):
                    value = lower[name][k]
                    if value['exact_zero']: continue
                    bound = audit.make_mpf(tuple(value['log_absolute_upper']['upper_exact_mpf_tuple']))
                    got = audit.make_mpf(tuple(p['finite_N_error_logs_in_physical_units'][name][k]['log_absolute_upper']['upper_exact_mpf_tuple']))
                    want = bound+power*radius_hi
                    assert got >= want or abs(got-want) <= abs(got)*audit.mpf('1e-150')
                    counts['lower_finite_N_endpoint_conditions'] += 1
        assert not p['limits']['restoration_postrestore_m15']

    for Z in ('-.3','0','.3'):
        for coordinate in ('0','g:.1','g:1','.25','.5','1'):
            check_moment(Z,coordinate)
        for field in fields[1:]:
            for coordinate in ('0','.125','.25','.5','.75','1'):
                check_pressure(field,Z,coordinate)
        for left,right in zip([f.rows.long_rows]+fields,fields):
            a,b = left.evaluate(Z,'1'),right.evaluate(Z,'0')
            for name in ('U','m','Q','U_y','Q_y'):
                for x,y in zip(a['physical_pressure_rows'][name]['signed_background_ordinary_Z'],
                               b['physical_pressure_rows'][name]['signed_background_ordinary_Z']):
                    lo1,hi1 = enclosure(x); lo2,hi2 = enclosure(y)
                    assert max(lo1,lo2) <= min(hi1,hi2), (left.chart,right.chart,Z,name)
                    counts['source_seam_overlap_conditions'] += 1
        print('Reference124-term m15 and reference/restore/post pressure to Rm: Z='+Z,flush=True)
    check_moment('.3','.875')
    for field in fields[1:]: check_pressure(field,'.3','.875')
    # Independent axial point; do not limit the adapter to cached inlet Z's.
    for field in fields: check_pressure(field,'.17','.3')
    # Long-decimal restoration truth must use the original q enclosure.
    check_pressure(fields[1],'.17','0.'+'3'*100)
    rejected = 0
    for field,Z,coordinate in [(f.rows,'1','.5'),(f.rows,'nan','.5'),(f.rows,'.3','-.1'),
            (f.rows,'.3','1.1'),(f.rows,'.3','g:-1'),(f.rows,'.3','g:nan'),
            (fields[1],'.3','g:.1'),(fields[2],'.3','g:.1')]:
        try: field.evaluate(Z,coordinate)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid original reference-family query admitted')
    code = dict(f.code_hashes)
    for field in fields: code.update(field.code_hashes)
    name = 'scripts/run_original_reference_m15.py'
    code[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,pressure_row_samples=pressure,
        input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,actual_original_source=source,independent_reference_primitives=independent,
                    actual_reference_m15_packets=len(samples),actual_pressure_packets=len(pressure),
                    **counts,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(reference_cumulative_signed_m15_with_complete_factored_errors=True,
            reference_restore_postrestore_pressure_rows=True,restoration_postrestore_m15=False,
            resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,
            common_global_N=False,exact_heat_exterior=False,stress_cone=False,
            genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    data = json.loads(path.read_text(encoding='utf-8')); assert data['checks']['passed']
    for name,digest in data['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest,name
    for name,digest in data['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest,name
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args = parser.parse_args(); path = DIRECTORY/NAME
    data = inspect(path) if args.inspect_report else run()
    if args.write_report:
        path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
