"""Focused buffer ODE/period proof checks and saved nonlinear integral replay.

Replays the saved actual inverse outputs, not expensive upstream inverses or
owners. Independently tests new reflection and weighted-bin/error lemmas.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_buffer_period_transport as source
import lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator as stress
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as compare
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv


def serialized_equal(actual,saved,path):
    exact_replay_equal(json.loads(json.dumps(source.encode(actual))),saved,path)


def independent_period_checks(family):
    c=MPIntervalContext();c.dps=100;p=mp.mp.clone();p.dps=350
    coordinates=source.native.HalfPstarCoordinates(c,c.mpf(0),family);t=coordinates.scalar(1)
    count=0;weight_count=0;error_count=0
    for A0 in ('-1.2','-.1','0','.1','1.2','1e-100'):
        for B0 in ('-.3','.2'):
            E=t.scalar(c.mpf('1.3'));A=t.scalar(c.mpf(A0));B=t.scalar(c.mpf(B0))
            got=source.reflected_density_mean(E,dict(A=A,B_over_Pstar=B),source.N)
            def density(sign):
                e=p.mpf('1.3');dE=e*p.expm1(sign*p.mpf(A0)/source.N);dV=sign*p.mpf(B0)/source.N
                return dict(m=dV,h=dE,k=e*dV+dE*dV,e=dV*dV-e*dE-dE*dE/2,p=e*dE+dE*dE/2)
            plus,minus=density(1),density(-1)
            for key in source.KEYS:compare.enclosed(p,c,got[key],(plus[key]+minus[key])/2);count+=1
    # Direct analytic weighted occupation masses, including canonical bins
    # wrapped by three arbitrary actual phase origins.
    N=16;w=p.mpf('.5');bins=8
    for rate in (p.mpf(0),p.mpf(1),p.mpf('1.5')):
        total=w if not rate else -p.expm1(-rate*w)/rate
        for origin in (p.mpf(0),p.mpf('.173'),p.mpf('.8')):
            for index in range(bins):
                mass=p.mpf(0)
                for cycle in range(-1,int(N*w)+2):
                    a=max(p.mpf(0),(cycle+p.mpf(index)/bins-origin)/N)
                    b=min(w,(cycle+p.mpf(index+1)/bins-origin)/N)
                    if a<b:mass+=(b-a) if not rate else (p.exp(-rate*(w-b))-p.exp(-rate*(w-a)))/rate
                assert total/bins*p.exp(-rate/N)<=mass<=total/bins*p.exp(rate/N);weight_count+=1
    # A varying-envelope periodic density has an explicit frozen mean.
    # Compare its actual weighted integral with that mean plus the proposed
    # full fixed-phase slow derivative and twice-oscillation error.
    p.dps=100
    N=8;w=p.mpf('.5');s=p.mpf('.2');Delta=p.mpf(1)/N
    for rate in (p.mpf(0),p.mpf(1),p.mpf('1.5')):
        total=w if not rate else -p.expm1(-rate*w)/rate
        for origin in (p.mpf('.17'),p.mpf('.61')):
            K=lambda y:p.exp(-y/2)*(s*p.sin(2*p.pi*(origin+N*y))+s*s*p.cos(2*p.pi*(origin+N*y))**2)
            points=[j*Delta for j in range(int(N*w)+1)]
            actual=p.quad(lambda y:p.exp(-rate*(w-y))*K(y),points)
            frozen=p.mpf(0)
            for j in range(int(N*w)):
                a=j*Delta;b=(j+1)*Delta
                mass=(b-a) if not rate else (p.exp(-rate*(w-b))-p.exp(-rate*(w-a)))/rate
                frozen+=mass*p.exp(-a/2)*s*s/2
            Kmax=s+s*s;Ly=Kmax/2
            error=total*(Delta*Ly+2*Kmax*p.expm1(rate*Delta))
            assert abs(actual-frozen)<=error;error_count+=1
    return dict(passed=True,independent_original_nonlinear_reflection_pair_comparisons=count,
        independent_arbitrary_origin_wrapped_weighted_phase_bin_checks=weight_count,
        independent_varying_envelope_and_period_weight_error_checks=error_count,
        tiny_nonzero_cosh_minus_one_factor_checked=True,finite_fixtures_not_actual_source_values=True)


def independent_buffer_ODE_checks(family):
    c=MPIntervalContext();c.dps=100;direct=MPIntervalContext();direct.dps=160;p=mp.mp.clone();p.dps=250
    coordinates=source.native.HalfPstarCoordinates(c,2*(c.exp(1)+11),family);template=coordinates.scalar(1)
    Z=['.37','.37'];z=source.Taylor(template,[c.mpf(Z),1,0,0,0,0]);qi=(1+z*z).reciprocal()
    EY=qi*source.prior.ScaledEnclosure(source.prior.FormalScale(template.scale.bases,(0,0,0,-1,0)),c.exp(c.mpf('5.8')),template.ledger)
    initial=dict(m=1+z,h=qi*c.mpf('.03'),k=qi*z*c.mpf('.02'),e=z*z*c.mpf('.001'),p=qi*qi*c.mpf('.1'))
    P0=z*z*c.mpf('.03')+c.mpf('.2')
    data=dict(exact_Z_range=Z,original_logC=source.encode(c.mpf('.1')),original_dstar_log=source.encode(c.mpf('-1')),
        original_parameter_sources=dict(original_eta_log=source.encode(c.mpf('-50'))),
        original_endpoint_parent_binding=dict(original_saved_native_inputs=dict(original_delta=source.encode(c.mpf('.001'))),
            original_separate_P0=source.native.records(P0),original_Utheta_endpoint_values=dict(right=source.native.records(EY)),
            original_right_history_values=source.native.records(initial)))
    adapter=source.OriginalBuffer(coordinates,data);count=0;root_checks=0
    zd=IntervalTaylor.variable(direct,direct.mpf(Z),5);qd=(1+zd*zd).reciprocal();logP=direct.exp(1)+11;P=direct.exp(logP)
    Eyd=qd*direct.exp(direct.mpf('5.8')-logP/2)
    old=dict(m=1+zd,h=qd*direct.mpf('.03'),k=qd*zd*direct.mpf('.02'),e=zd*zd*direct.mpf('.001'),p=qd*qd*direct.mpf('.1'))
    P0d=zd*zd*direct.mpf('.03')+direct.mpf('.2');zero=zd*0
    for a,b in (('0','.5'),('5','5.5'),('10.5','11')):
        box=c.mpf((a,b));got=adapter.source(box);roots=adapter.signed_roots(got,box)['roots']
        for s in (direct.mpf(a),(direct.mpf(a)+direct.mpf(b))/2,direct.mpf(b)):
            dec=direct.exp(-s);d3=direct.exp(-direct.mpf('1.5')*s);E=Eyd*direct.exp(-s/2)
            hist=dict(m=old['m']*dec,h=old['h']*d3+Eyd*(direct.exp(-s/2)-d3),k=old['k']*d3,
                e=(old['e']-Eyd*Eyd*(s/2))*dec,p=old['p']+Eyd*Eyd*((1-dec)/2))
            physical=source.native.original.physical_mixed(direct,Z,direct.mpf('.001'),E,[zero-direct.mpf('.5')]+[zero]*3,[zero]*5,hist,P0d,1/P**2)
            for name in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4'):
                for key,rows in physical[name].items():
                    for order,value in rows.items():
                        for x in ep(value):compare.enclosed(p,c,got['physical'][name][key][order],p.mpf(x));count+=1
            U=physical['physical_velocity_pressure_y_derivative_Taylor']['Utheta_over_Pstar']
            H=physical['actual_normalized_primitive_y_derivative_axial5'];press=[P0d+H['p'][0]]+H['p'][1:]
            parts=stress.raw_pre_stress_rows(direct,direct.mpf('.001'),zd,U,[zero]*5,H,press)
            radius=direct.exp(10*direct.mpf('.1')+direct.ln(110)+11*logP-11+s)
            for label,omit in (('theta','variable_radial_shear'),('axial','axial_radial_shear')):
                rows=[zero]*3
                for name,part in parts[label].items():
                    if name==omit:continue
                    factor=1 if label=='theta' else 1/P if part['mode'][1]==0 else P
                    rows=[rows[j]+part['shape'][j]*factor for j in range(3)]
                nums=[sum((rows[i]*math.comb(j,i) for i in range(j+1)),zero)*radius for j in range(3)]
                ratios=[nums[0]/U[0]]
                ratios.append((nums[1]-ratios[0]*U[1])/U[0])
                ratios.append((nums[2]-ratios[0]*U[2]-ratios[1]*U[1]*2)/U[0])
                for j in range(3):
                    for n in (0,1):
                        for x in ep(ratios[j][n]):compare.enclosed(p,c,roots['p1' if label=='theta' else 'p2'][(j,n)],p.mpf(x));root_checks+=1
    return dict(passed=True,independent_original_buffer_physical_primitive_mixed4_comparisons=count,
        independent_original_signed_stress_ordinary_y2_Z1_unit_radius_comparisons=root_checks,
        original_ordinary_ODE_and_stress_programs_independently_evaluated=True,finite_M1_fixture_not_actual_M40_source=True)


def replay_saved_cell(adapter,row):
    c=adapter.c;coordinates=adapter.coordinates;rec=row['actual_period_source_and_integrals']
    endpoints=[Fraction(x) for x in row['original_exact_t_endpoints']]
    cv=lambda f:c.mpf(f.numerator)/f.denominator
    a,b=map(cv,endpoints);box=c.mpf((a,b));width=b-a;bins=source.PHASE_BINS
    got=adapter.source(box);roots=adapter.signed_roots(got,box)
    serialized_equal(got['record'],rec['original_source'],'buffer-source')
    serialized_equal(roots['record'],rec['original_signed_stress_source'],'buffer-stress')
    serialized_equal(source.phase_origin(adapter),rec['actual_original_phase_binding'],'actual-origin')
    restore=lambda v:source.complete.restore_half_source(v,coordinates)
    bybin={i:[] for i in range(bins)};phase_queries=0
    for query in rec['actual_original_nonlinear_phase_density_queries']:
        index=query['original_phase_bin'];inverse=query['original_conditional_first_jets']
        assert inverse['status']=='enclosed' and inverse['slow_derivatives_hold_actual_phi_fixed']
        phi=c.mpf((c.mpf(index)/bins,c.mpf(index+1)/bins))
        assert ep(iv(c,query['actual_complete_period_phase_box']))==ep(phi)
        values={k:restore(v) for k,v in inverse['original_A_B_first_derivative_enclosures'].items()}
        E=roots['roots']['E'][source.ZERO];EZ=roots['roots']['E'][(0,1)];Ey=roots['roots']['E'][(1,0)];zero=E.scalar(0)
        den=adapter.kernel(E,EZ,zero,zero,values,source.N)
        slow=adapter.kernel(E,Ey,zero,zero,dict(values,A_Z=values['A_y'],B_Z_over_Pstar=values['B_y_over_Pstar']),source.N)['Z_derivatives']
        mean=source.reflected_density_mean(E,values,source.N) if index<bins//2 else None
        serialized_equal(source.native.records(den['kernels']),query['actual_nonlinear_five_C0_density'],'saved-C0-density')
        serialized_equal(source.native.records(den['Z_derivatives']),query['actual_nonlinear_five_Z_density'],'saved-Z-density')
        serialized_equal(source.native.records(slow),query['actual_complete_fixed_phase_slow_y_density'],'saved-slow-y-density')
        serialized_equal(None if mean is None else source.native.records(mean),query['actual_reflected_C0_pair_density'],'saved-reflection-density')
        bybin[index].append(dict(density=den,slow=slow,mean=mean));phase_queries+=1
    assert all(bybin.values())
    union=source.rc.density.local.same_source_union;zero=coordinates.scalar(0)
    density=[{k:union([d['density']['kernels'][k] for d in bybin[i]]) for k in source.KEYS} for i in range(bins)]
    slow=[{k:union([d['slow'][k] for d in bybin[i]]) for k in source.KEYS} for i in range(bins)]
    jets=[{k:union([d['density']['Z_derivatives'][k] for d in bybin[i]]) for k in source.KEYS} for i in range(bins)]
    means=[{k:union([d['mean'][k] for d in bybin[i]]) for k in source.KEYS} for i in range(bins//2)]
    whole={k:union([d[k] for d in density]) for k in source.KEYS};whole_y={k:union([d[k] for d in slow]) for k in source.KEYS}
    avg={k:sum((d[k]*(c.mpf(2)/bins) for d in means),zero) for k in source.KEYS}
    for value,key in ((whole,'complete_fixed_phase_density_C0_cover'),(whole_y,'complete_fixed_phase_slow_y_density_cover'),(avg,'fixed_slow_reflection_mean_C0_cover')):
        serialized_equal(source.native.records(value),rec[key],key)
    geo=dict(width=coordinates.scalar(width),regular=width,scalar_cover=width,record=rec['original_true_geometry'])
    values={};Zvalues={}
    for key,rate in source.rc.RATES.items():
        r=Fraction(rate);rr=c.mpf(r.numerator)/r.denominator
        factor=source.rc.transfer.true_width_kernel(coordinates,geo,rate)
        saved=rec['original_true_width_Duhamel_factors'][key]
        serialized_equal(dict(mass=factor['mass'].record(),decay=factor['decay'].record(),branch=factor['branch']),saved,'Duhamel-factor')
        slow_error=source.absolute_cover(whole_y[key])*(c.mpf(1)/source.N)
        weight_error=source.absolute_cover(whole[key])*(2*c.expm1(rr/source.N))
        values[key]=factor['mass']*(avg[key]+source.symmetric_cover(slow_error+weight_error))
        weights=c.exp(c.mpf((-ep(rr/source.N)[1],ep(rr/source.N)[1])))*(c.mpf(1)/bins)
        Zvalues[key]=factor['mass']*sum((d[key]*weights for d in jets),zero)
        serialized_equal(dict(full_slow_envelope_error=slow_error.record(),period_weight_error=weight_error.record()),rec['original_per_rate_period_and_slow_errors'][key],'explicit-errors')
    serialized_equal(source.native.records(values),rec['actual_local_five_C0_integrals'],'C0-integral')
    serialized_equal(source.native.records(Zvalues),rec['actual_local_five_Z_integrals'],'Z-integral')
    assert rec['primitive_support_caps_not_used_as_field_values'] and rec['C0_reflection_mean_is_frozen_slow_with_explicit_errors']
    assert rec['Z_direct_density_phase_integrals_not_exact_reflection_cancellation']
    return dict(geometry=geo,values=values,Z_derivatives=Zvalues),phase_queries


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N and manifest['exact_periods_per_buffer']==11*source.N
    assert manifest['original_whole_buffer_C0_Z_integrals_and_new_Rd_Rc_histories_executed']
    assert not manifest['actual_functional_targets_controls_full_Z_axis_N_heat_stress_recursion_admitted']
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    assert all(manifest[k] is False for k in source.common.current.FLAGS)
    c=MPIntervalContext();c.dps=240;cells=0;phase_queries=0
    with mp.workdps(300):
        independent=independent_period_checks(manifest['source_family']);ode=independent_buffer_ODE_checks(manifest['source_family'])
        for archive in manifest['actual_original_buffer_period_transport_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            data=json.loads(raw);parents=[]
            for key in ('accepted_original_axial_archive','accepted_original_transition_archive'):
                parent=data[key];praw=gzip.decompress((HERE/parent['filename']).read_bytes())
                assert hashlib.sha256(praw).hexdigest()==parent['lossless_original_json_sha256'];parents.append(json.loads(praw))
            original,transition=parents
            coordinates=source.native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),manifest['source_family'])
            adapter=source.OriginalBuffer(coordinates,original);operator=source.rc.history.C1DuhamelOperator(coordinates)
            rows=data['actual_original_buffer_source_rows'];assert len(rows)==source.CELL_COUNT
            for index,row in enumerate(rows):
                assert list(map(Fraction,row['original_exact_t_endpoints']))==[Fraction(11*index,source.CELL_COUNT),Fraction(11*(index+1),source.CELL_COUNT)]
                local,n=replay_saved_cell(adapter,row);phase_queries+=n;cells+=1
                source.rc.transfer.append_true_cell(operator,local['geometry'],local['values'],local['Z_derivatives'],manifest['source_family'])
                print('Saved original buffer density/integral replay:',data['exact_Z_range'],index+1,flush=True)
            operator.coefficients={k:coordinates.decay(11,r) for k,r in source.rc.RATES.items()}
            update=source.propagate_to_Rc(coordinates,original,operator,transition)
            serialized_equal(update,data['actual_updated_Rd_Rc_correction_and_own_history'],'genuine-new-Rd-Rc')
            assert data['complete_original_buffer_C0_Z_source_and_operator_executed']
            assert not data['full_original_target_controls_N_axis_heat_stress_recursion_admitted']
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_nonlinear_period_averaging_checks=independent,independent_original_buffer_ODE_and_stress_checks=ode,
        actual_complete_buffer_source_cells_replayed=cells,actual_saved_phase_inverse_nonlinear_density_queries_replayed=phase_queries,
        actual_five_C0_Z_local_integral_rows_replayed=cells*10,exact_11_unit_coefficients_and_genuine_new_Rd_Rc_histories_replayed=True,
        nonlinear_density_recomputed_from_lossless_original_inverse_outputs=True,expensive_actual_inverse_solvers_not_repeated=True,
        accepted_upstream_or_transition_producers_not_rerun=True,
        frozen_C0_mean_has_explicit_full_slow_and_weight_errors_Z_not_claimed_exact_cancellation=True,
        actual_five_functional_defects_controls_full_Z_axis_N_heat_stress_recursion_still_open=True,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Complete original buffer period integration and new Rd/Rc C1 histories PASS; functional closure open',flush=True)
    return result


if __name__=='__main__':run()
