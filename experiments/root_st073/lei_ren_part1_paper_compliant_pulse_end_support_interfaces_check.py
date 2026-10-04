"""Focused four-support full-stress/error interfaces and retained-history fixture."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_end_support_interfaces import (
    source_difference_proof, full_difference_rows, difference_transport,
    symmetric_jet, beta_tail_bound, EDGES, PREFIX, FALSE_FLAGS)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import pulse_coefficients
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import pulse_velocity_rows,pulse_remainder_sectors
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent


def retained_history_fixture():
    """Independent original-beta integrals on either side of each exact edge."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=100
        ell=mp.mpf('.15');mu=mp.mpf('.1');delta=mp.mpf('.06');Z=mp.mpf('.31')
        raw=lambda r:mp.exp(-1/(1-r*r)) if abs(r)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1]);cj=IntervalTaylor(c,['-.2','.05','.02',0,0,0])
        z=IntervalTaylor.variable(c,mp.nstr(Z,90),5)
        # Polynomial control is written in powers of Z-Z0; whole axial5 jets retained.
        C=IntervalTaylor(c,[1+c.mpf(mp.nstr(Z,90))**2,2*c.mpf(mp.nstr(Z,90)),1,0,0,0]).reciprocal()
        M0=IntervalTaylor(c,['.7','.1','.02',0,0,0])
        N0=IntervalTaylor(c,['.4','-.03','.01',0,0,0])
        J0=IntervalTaylor(c,['.2','.04','.01',0,0,0])
        zero=C*0;counts={};worst=mp.mpf(0);tol=mp.mpf('1e-55')
        def point(value):return c.mpf(mp.nstr(value,90))
        def midpoint(value):
            lo,hi=endpoints(value);return (lo+hi)/2
        def compare(label,bound,value):
            nonlocal worst
            lo,hi=endpoints(bound);miss=max(lo-value,value-hi,mp.mpf(0));worst=max(worst,miss)
            if miss>tol*max(1,abs(value)):raise ArithmeticError('Retained-history flat difference failed '+label)
            counts[label.split('/')[0]]=counts.get(label.split('/')[0],0)+1
        nonzero_histories=0
        for edge in EDGES:
            edgevalue=mp.mpf(edge['center'])+edge['side']*ell
            for h in (mp.mpf('.01'),mp.mpf('.003')):
                sv=edgevalue-edge['side']*h
                beta=lambda value:raw((value-edge['center'])/ell)/(ell*norm)
                B=[cj*point(mp.diff(beta,sv,k)) for k in range(5)]
                def integral(rate,power):
                    if sv>=edgevalue:
                        return mp.quad(lambda t:mp.exp(-rate*(sv-t))*beta(t)**power,[edgevalue,(edgevalue+sv)/2,sv])
                    return -mp.quad(lambda t:mp.exp(-rate*(sv-t))*beta(t)**power,[sv,(edgevalue+sv)/2,edgevalue])
                dm=[cj*point(integral(mp.mpf('.5')-mu,1))]
                dn=[cj*point(integral(mp.mpf('.5')-2*mu,1))]
                dj=[-(cj*cj)*point(integral(-2*mu,2))]
                for k in range(4):
                    dm.append(B[k]-dm[k]*point(mp.mpf('.5')-mu))
                    dn.append(B[k]-dn[k]*point(mp.mpf('.5')-2*mu))
                    square=sum((B[l]*B[k-l]*math.comb(k,l) for l in range(k+1)),zero)
                    dj.append(dj[k]*point(2*mu)-square)
                rm=[M0*point(mp.exp(-(mp.mpf('.5')-mu)*(sv-edgevalue)))]
                rn=[N0*point(mp.exp(-(mp.mpf('.5')-2*mu)*(sv-edgevalue)))]
                rj=[J0*point(mp.exp(2*mu*(sv-edgevalue)))]
                for k in range(4):
                    rm.append(rm[k]*point(-(mp.mpf('.5')-mu)))
                    rn.append(rn[k]*point(-(mp.mpf('.5')-2*mu)))
                    rj.append(rj[k]*point(2*mu))
                am=[rm[k]+dm[k] for k in range(5)]
                an=[rn[k]+dn[k] for k in range(5)]
                aj=[rj[k]+dj[k] for k in range(5)]
                e=[C*0+mp.mpf('.8')];p=[-C*C*point(mp.mpf('.4'))]
                for k in range(4):
                    e.append(e[k]*point(2*mu)-(zero+mp.mpf('.5') if k==0 else zero))
                    p.append(p[k]*point(1+2*mu)+(C*C/2 if k==0 else zero))
                params=(point(delta),point(mu),z,C,point(mp.mpf('.3')))
                actual=pulse_coefficients(*params,B,am,an,e,aj,p)
                reference=pulse_coefficients(*params,[zero]*5,rm,rn,e,rj,p)
                forcing=[symmetric_jet(cj*(beta_tail_bound(c,k,2*point(h)/point(ell))/(point(ell)**(k+1)*point(norm)))) for k in range(5)]
                history=difference_transport(c,forcing,point(mu),point(h))
                bounds=full_difference_rows(c,point(delta),point(mu),z,C,point(mp.mpf('.3')),forcing,
                    history['linear_m1'],history['linear_m2'],[-v for v in history['energy']],am)
                for label in ('theta','axial'):
                    for name,sector in actual[label].items():
                        for k in range(4):
                            for n in range(4-k):
                                target=midpoint(sector['full_derivative_rows'][k][n])-midpoint(reference[label][name]['full_derivative_rows'][k][n])
                                compare('stress/'+label+'/'+name+'/'+str(k)+str(n),bounds['stress'][label][name]['full_derivative_rows'][k][n],target)
                av=pulse_velocity_rows(c,point(delta),point(mu),z,C,B,am)
                rv=pulse_velocity_rows(c,point(delta),point(mu),z,C,[zero]*5,rm)
                ae=pulse_remainder_sectors(c,point(delta),point(mu),z,av)
                re=pulse_remainder_sectors(c,point(delta),point(mu),z,rv)
                for label in ('radial','theta','axial'):
                    for k in range(5):
                        for n in range(5-k):
                            compare('velocity/'+label+'/'+str(k)+str(n),bounds['velocity'][label][k][n],midpoint(av[label][k][n])-midpoint(rv[label][k][n]))
                    for name,sector in ae[label].items():
                        for k in range(3):
                            for n in range(3-k):
                                compare('error/'+label+'/'+name+'/'+str(k)+str(n),bounds['remainder'][label][name]['rows'][k][n],
                                    midpoint(sector['rows'][k][n])-midpoint(re[label][name]['rows'][k][n]))
                if abs(midpoint(rv['radial'][0][0]))<mp.mpf('.01'):raise ArithmeticError('Fixture lost nonzero reference Ur')
                nonzero_histories+=1
        return dict(all_passed=True,checks=counts,total_checks=sum(counts.values()),tolerance=str(tol),
            maximum_positive_enclosure_miss=mp.nstr(worst,25),
            independent_original_beta_quadrature_and_derivatives=True,both_edge_orientations_and_all_four_edges=True,
            actual_nonzero_reference_radial_histories=nonzero_histories,
            fixture_parameters_not_actual_source_values=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'pulse_end_support_interfaces.json';record=json.loads((HERE/name).read_bytes())
        hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Interface source changed '+path)
        if record['source_difference_proof']!=source_difference_proof():raise ValueError('Interface proof changed')
        c=MPIntervalContext();c.dps=300;finite=zeros=0;groups={}
        def bounds(point):
            values=[]
            for jet in point['unscaled_beta_forcing_difference_rows']:values.extend(jet['coefficients'])
            for jets in point['source_primitive_differences'].values():
                for jet in jets:values.extend(jet['coefficients'])
            for name in ('similarity_stress_difference_mixed3','physical_stress_difference_mixed3_coefficients',
                'physical_three_component_error_difference_mixed2_coefficients'):
                for sectors in point[name].values():
                    for grid in sectors.values():values.extend(grid.values())
            for grid in point['similarity_velocity_difference_mixed4'].values():values.extend(grid.values())
            return values
        for point in record['interfaces']:
            key=(point['edge']['row'],point['edge']['side']);groups.setdefault(key,[]).append(point)
            atzero=endpoints(read_interval(c,point['h']))==(0,0)
            for value in bounds(point):
                lo,hi=endpoints(read_interval(c,value))
                if not all(mp.isfinite(v) for v in (lo,hi)):raise ArithmeticError('Nonfinite flat difference')
                finite+=1
                if atzero:
                    if lo or hi:raise ArithmeticError('Zero distance difference lost exact zero')
                    zeros+=1
            for flag in ('actual_boundary_histories_retained','exact_shared_swirl_pressure_and_incoming_memory_differences_zero',
                'current_positive_log_factors_preserved_as_parent_recipes'):
                if not point[flag]:raise ValueError('Current history/source omitted')
        if len(groups)!=4 or any(len(points)!=4 for points in groups.values()):raise ValueError('Original four support edges not fully covered')
        monotone=0
        for points in groups.values():
            norms=[[max(abs(v) for v in endpoints(read_interval(c,value))) for value in bounds(point)] for point in points]
            for left,right in zip(norms,norms[1:]):
                for a,b in zip(left,right):
                    if b>a:raise ArithmeticError('Flat interface envelope grew toward endpoint')
                    monotone+=1
        for flag in FALSE_FLAGS:
            if record[flag]:raise ValueError('Interface scope overclaimed '+flag)
        current=json.loads((HERE/(PREFIX+'pulse_end_physical_C2_check.json')).read_bytes())
        if not current['all_passed'] or not current['independent_full_meridional_Cartesian_oracle']['all_passed']:
            raise ValueError('Accepted full physical operator required')
        print('Four whole-Z support boundaries: exact-zero, finite and vanishing mixed bounds PASS',flush=True)
        fixture=retained_history_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
            implicit_source_sha256=record['implicit_source_sha256'],
            actual_pulse_end_all_four_support_functional_interfaces_verified=True,
            actual_pulse_end_stress3_and_physical_error2_flat_interface_bounds_available=True,
            finite_difference_bounds_checked=finite,exact_endpoint_difference_zeros=zeros,monotone_bound_comparisons=monotone,
            source_interface_identities=len(record['source_difference_proof']['identities']),
            retained_nonzero_history_fixture=fixture,
            accepted_current_full_Cartesian_operator_consumed_without_rerun=True,
            full_stress3_implies_completed_diagonal2_divergence2_flat_joins=True,
            support_flatness_does_not_prove_global_temporal_flat_remainder=True,
            **{flag:False for flag in FALSE_FLAGS},input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS four original support full-stress/physical-error interfaces with retained nonzero histories; cone/global pending',flush=True)
        return result


if __name__=='__main__':run()
