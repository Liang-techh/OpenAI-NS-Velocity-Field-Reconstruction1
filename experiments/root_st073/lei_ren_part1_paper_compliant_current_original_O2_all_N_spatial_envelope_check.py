"""Independent unpaired remainder/IBP diagnostics and actual source check."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_all_N_spatial_envelope as current
import lei_ren_part1_paper_compliant_current_original_O2_uniform_mean_bias_check as previous

saved=previous.saved;ep=current.ep;C0,Y,Z,YZ=current.mixed.ORDERS


def independent_full_unpaired_remainder(owner):
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas;c=owner.ctx
    p=mp.mp.clone();p.dps=180;z0=p.mpf('-.37');comparisons=0;cases=0
    for N in (160,257,10**9):
        for vv in ('-.21','.21'):
            for aa in ('-.5','0','1e-50','.5'):
                data=dict(E=('1.1','-.3'),V=(vv,'.4'),A=(aa,'-.13'),B=('.53','.23'))
                jets={name:current.MixedJet(a,{order:a.scalar(c.mpf(row[0] if order==C0 else row[1] if order==Z else 0))
                    for order in current.mixed.ORDERS}) for name,row in data.items()}
                def f(name,z):
                    value,derivative=map(p.mpf,data[name]);return value+derivative*(z-z0)
                def remainder(key,z):
                    E,V,A,B=(f(name,z) for name in ('E','V','A','B'))
                    de=E*p.expm1(A/N);dv=B/N
                    delta=dict(m=dv,h=de,k=V*de+E*dv+de*dv,
                        e=2*V*dv+dv*dv-E*de-de*de/2,p=E*de+de*de/2)
                    first=dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)
                    return N*N*(delta[key]-first[key]/N)
                targets=[{key:remainder(key,z0) for key in current.integrals.KEYS},
                    {key:p.diff(lambda z:remainder(key,z),z0) for key in current.integrals.KEYS}]
                for physical in (False,True):
                    caps=current.mean.native_coefficients(a,**jets,physical_amplitude=physical)
                    upper=lambda group,key:ep(caps[group][key].finite_interval())[1]
                    values=dict(m=c.mpf(0),h=c.mpf((0,upper('C0','h'))),k=c.mpf((-upper('C0','k'),upper('C0','k'))),
                        e=c.mpf((-upper('C0','e_negative'),upper('C0','e_positive'))),p=c.mpf((0,upper('C0','p'))))
                    derivatives={key:c.mpf((-upper('Z',key),upper('Z',key))) for key in current.integrals.KEYS}
                    for group,targetsrow in zip((values,derivatives),targets,strict=True):
                        for key,target in targetsrow.items():saved.contains(group[key],target,p.mpf('1e-130'));comparisons+=1
                cases+=1
    assert current.exact_spatial_theorem()['exact_full_changed_density_first_order_second_remainder_identities']==5
    return dict(passed=True,compatible_unpaired_remainder_cases=cases,
        independent_full_exponential_remainder_C0_Z_comparisons=comparisons,
        N_values=[160,257,10**9],both_V_signs_and_full_phase_A_boundary_zero_tiny=True,
        mean_bias_retained_in_unpaired_remainder=True)


def independent_nonzero_endpoint_spatial_IBP(owner):
    """Compatible analytic fields, true rapid phase and independent quadrature."""
    p=mp.mp.clone();p.dps=65;c=owner.ctx;left,right=p.mpf(12)/64,p.mpf(13)/64;z0=p.mpf('.17')
    identitychecks=0;integralchecks=0;endpointchecks=0
    def parameters(y,z):
        E=1+p.mpf('.1')+p.mpf('.02')*y-p.mpf('.03')*z+p.mpf('.01')*y*z
        V=p.mpf('.21')+p.mpf('.04')*y+p.mpf('.02')*z
        A=p.mpf('.3')+p.mpf('.05')*y+p.mpf('.01')*z+p.mpf('.02')*y*z
        B=p.mpf('.53')-p.mpf('.04')*y+p.mpf('.03')*z-p.mpf('.01')*y*z
        return E,V,A,B
    def amplitude(key,y,z):
        E,V,A,B=parameters(y,z)
        return dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)[key]
    # Honest polynomial source and derivative covers over one finite cell.
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas
    yy=c.mpf((str(left),str(right)));zz=c.mpf(str(z0));sine=c.mpf((-1,1))
    rows=dict(E=(1+c.mpf('.1')+c.mpf('.02')*yy-c.mpf('.03')*zz+c.mpf('.01')*yy*zz,c.mpf('.02')+c.mpf('.01')*zz,-c.mpf('.03')+c.mpf('.01')*yy,c.mpf('.01')),
        V=(c.mpf('.21')+c.mpf('.04')*yy+c.mpf('.02')*zz,c.mpf('.04'),c.mpf('.02'),c.mpf(0)),
        A=((c.mpf('.3')+c.mpf('.05')*yy+c.mpf('.01')*zz+c.mpf('.02')*yy*zz)*sine,(c.mpf('.05')+c.mpf('.02')*zz)*sine,(c.mpf('.01')+c.mpf('.02')*yy)*sine,c.mpf('.02')*sine),
        B=((c.mpf('.53')-c.mpf('.04')*yy+c.mpf('.03')*zz-c.mpf('.01')*yy*zz)*sine,(-c.mpf('.04')-c.mpf('.01')*zz)*sine,(c.mpf('.03')-c.mpf('.01')*yy)*sine,-c.mpf('.01')*sine))
    jets={name:current.MixedJet(a,{order:a.scalar(value) for order,value in zip(current.mixed.ORDERS,row,strict=True)}) for name,row in rows.items()}
    leading=current.slow.leading_densities(a,**jets);rc=current.mean.native_coefficients(a,**jets,physical_amplitude=True)
    for N in (160,257):
        phase0=p.mpf('.137');points=[left,right]
        start=int(p.floor(N*left+phase0));end=int(p.ceil(N*right+phase0))
        for k in range(start,end+1):
            x=(k-phase0)/N
            if left<x<right:points.append(x)
        points=sorted(points)
        def G(key,y,z):return amplitude(key,y,z)*(1-p.cos(2*p.pi*(N*y+phase0)))/(2*p.pi)
        def Gy(key,y,z):return p.diff(lambda v:amplitude(key,v,z),y)*(1-p.cos(2*p.pi*(N*y+phase0)))/(2*p.pi)
        def full(key,y,z):
            E,V,A,B=parameters(y,z);sin=p.sin(2*p.pi*(N*y+phase0));de=E*p.expm1(A*sin/N);dv=B*sin/N
            return dict(m=dv,h=de,k=V*de+E*dv+de*dv,e=2*V*dv+dv*dv-E*de-de*de/2,p=E*de+de*de/2)[key]
        for key,rate in current.integrals.five.RATES.items():
            rr=p.mpf(rate);K=lambda y:p.exp(-rr*(1-y))
            mass=(right-left) if not rr else (K(right)-K(left))/rr
            for j,(order,dy) in enumerate(((C0,Y),(Z,YZ))):
                diff=lambda fun,y:fun(key,y,z0) if j==0 else p.diff(lambda z:fun(key,y,z),z0)
                endpoint=K(right)*diff(G,right)-K(left)*diff(G,left)
                transport=p.quad(lambda y:K(y)*(diff(Gy,y)+rr*diff(G,y)),points)
                direct=p.quad(lambda y:K(y)*(amplitude(key,y,z0) if j==0 else p.diff(lambda z:amplitude(key,y,z),z0))*p.sin(2*p.pi*(N*y+phase0))/N,points)
                assert abs(direct-(endpoint-transport)/N**2)<p.mpf('1e-55');identitychecks+=1
                assert abs(endpoint)>p.mpf('1e-20');endpointchecks+=1
                cap=lambda value:max(abs(v) for v in ep(value.finite_interval()))
                gcap=cap(leading[key][order])/2;gycap=cap(leading[key][dy])/2
                if j==0:
                    rest=rc['C0'][key] if key!='e' else a.add(rc['C0']['e_negative'],rc['C0']['e_positive'])
                else:rest=rc['Z'][key]
                bound=(K(left)+K(right))*gcap+mass*(gycap+rr*gcap+cap(rest))
                actual=p.quad(lambda y:K(y)*diff(full,y),points)
                assert N*N*abs(actual)<=bound+p.mpf('1e-55');integralchecks+=1
    return dict(passed=True,independent_true_spatial_C0_Z_IBP_identities=identitychecks,
        independent_full_nonlinear_spatial_contribution_bounds=integralchecks,nonzero_kept_endpoint_terms=endpointchecks,
        actual_rapid_phase_not_uniform_phase_average=True,slow_y_and_phase_transport_separated=True,
        finite_compatible_fields_not_original_source_substitution=True)


def actual_source_and_mass_reconstruction(owner,manifest):
    c=owner.ctx;p=mp.mp.clone();p.dps=150;live=owner.integrate();encoded=current.base.encoded(live['report'])
    report=manifest['actual_original_O2_full_spatial_all_N_envelope'];rows=report['actual_source_spatial_IBP_records']
    assert report['exact_y_window']==['0','1'] and report['exact_Z_range']==['-1','1']
    assert report['all_integer_N_lower']==160 and report['source_family']==owner.family and len(rows)==64
    assert report['own_rates']==current.integrals.five.RATES
    counts=dict(primitive_rows=0,remainder_rows=0,components=0,masses=0,kernel_endpoints=0,unions=0)
    for index,(row,new) in enumerate(zip(rows,encoded['actual_source_spatial_IBP_records'],strict=True)):
        assert row['exact_y_cell']==new['exact_y_cell'] and row['branch_overlap_hulled_not_added']
        assert row['all_each_cell_and_global_endpoints_retained'] and row['no_extra_R_Jacobian_period_count_or_second_mean_bias']
        assert {b['branch'] for b in row['conditional_original_source_spatial_IBP_records']}==set(current.integrals.BRANCHES)
        left,right=p.mpf(index)/64,p.mpf(index+1)/64
        for branch,other in zip(row['conditional_original_source_spatial_IBP_records'],new['conditional_original_source_spatial_IBP_records'],strict=True):
            frame=owner.source.frame(64,index,branch=branch['branch'])
            assert branch['actual_original_predicate']==frame.record['actual_original_source_predicate']
            primitive=branch['periodic_primitive'];assert primitive['actual_primitive_point_value_not_evaluated']
            assert primitive['phase_mean_zero_from_same_original_reflection'] and primitive['actual_slow_function_derivatives_used_not_cap_derivatives']
            assert mp.mp.make_mpf(tuple(primitive['positive_integration_length_bound']['exact_mpf_tuple']))==.5
            for key,orders in primitive['native_periodic_primitive_absolute_caps'].items():
                for order,value in orders.items():previous.exact_native(c,value,other['periodic_primitive']['native_periodic_primitive_absolute_caps'][key][order]);counts['primitive_rows']+=1
            remainder=branch['full_nonlinear_remainder'];assert remainder['mean_bias_not_assumed_zero_or_subtracted']
            assert remainder['uniform_exp_floor_uses_global_reflection_extended_A_bound']
            for family,again in zip(remainder['native_positive_all_N_remainder_families'],other['full_nonlinear_remainder']['native_positive_all_N_remainder_families'],strict=True):
                for group,othergroup in zip(family,again,strict=True):
                    for key,value in group.items():previous.exact_native(c,value,othergroup[key]);counts['remainder_rows']+=1
            for j,(group,again) in enumerate(zip(branch['original_C0_Z_IBP_components'],other['original_C0_Z_IBP_components'],strict=True)):
                for key,part in group.items():
                    for name in ('native_kept_both_endpoint_caps','native_slow_and_kernel_transport_cap','native_full_nonlinear_remainder_cap','native_total_Nminus2_coefficient'):
                        previous.exact_native(c,part[name],again[key][name]);counts['components']+=1
                    previous.exact_interval(c,part['normalized_dominating_unit_coefficient_cap'],again[key]['normalized_dominating_unit_coefficient_cap'])
                    if branch['branch']=='regular' and j==0:
                        rate=p.mpf(current.integrals.five.RATES[key]);ka=p.exp(-rate*(1-left));kb=p.exp(-rate*(1-right))
                        mass=right-left if not rate else (kb-ka)/rate
                        saved.contains(saved.interval(c,part['positive_own_rate_mass']),mass,p.mpf('1e-120'));counts['masses']+=1
                        saved.contains(saved.interval(c,part['left_kernel']),ka,p.mpf('1e-120'))
                        saved.contains(saved.interval(c,part['right_kernel']),kb,p.mpf('1e-120'));counts['kernel_endpoints']+=2
        for group,again in zip(row['normalized_C0_Z_Nminus2_coefficient_unions'],new['normalized_C0_Z_Nminus2_coefficient_unions'],strict=True):
            for key,value in group.items():previous.exact_interval(c,value,again[key]);counts['unions']+=1
    for group,again in zip(report['normalized_full_spatial_C0_Z_Nminus2_coefficient_caps'],encoded['normalized_full_spatial_C0_Z_Nminus2_coefficient_caps'],strict=True):
        for key,value in group.items():previous.exact_interval(c,value,again[key]);counts['unions']+=1
    for unit,newunit,order in zip(report['dominating_C0_Z_output_units'],encoded['dominating_C0_Z_output_units'],current.OUTPUT_ORDERS,strict=True):
        previous.exact_native(c,unit,newunit)
        assert tuple(unit['formal_positive_scale']['source_exponents'])+(unit['formal_positive_scale']['radius_power'],)==current.slow.UNITS[order]
    return live,dict(passed=True,**counts,complete_actual_source_cover=True,
        each_own_rate_mass_and_each_local_global_endpoint_preserved=True,
        actual_nonlinear_mean_bias_inside_remainder_not_added_twice=True)


def candidate_phase_units_and_guards(owner,live,manifest):
    c=owner.ctx;queries=0;phasechecks=0;rejected=0
    for text,stored in manifest['candidate_full_spatial_queries'].items():
        N=int(text);query=owner.at_candidate(live,N=N);new=current.base.encoded(query['record'])
        assert stored['candidate_not_global_N_admission'] and stored['approximate_phase_not_used_as_selected_value']
        for j,group in enumerate(stored['normalized_signed_C0_Z_contribution_envelopes']):
            for key,value in group.items():
                upper=ep(live['coefficients'][j][key]*(c.mpf(1)/(c.mpf(N)**2)))[1]
                previous.exact_interval(c,value,c.mpf((-upper,upper)));queries+=1
        for name in ('native_C0_contribution_envelopes','native_ordinary_Z_contribution_envelopes'):
            for key,value in stored[name].items():previous.exact_native(c,value,new[name][key]);queries+=1
        endpoints=stored['actual_original_radius_common_N_phase_endpoints'];assert len(endpoints)==65
        for index,(row,again) in enumerate(zip(endpoints,new['actual_original_radius_common_N_phase_endpoints'],strict=True)):
            assert row['original_y_exact']==str(current.s.Rational(index,64)) and row['explicit_candidate_N']==N
            assert row['exact_source_phase']=='frac(N*(log(110/4)+14*logPstar+10*logCstar+1000+y-hb*s_c/2))'
            assert row['source_family']==owner.family and row['candidate_N_not_global_selection']
            assert row['positive_original_origin_offset']['strictly_positive'] and row['positive_original_origin_offset']['actual_hb_sc_product_not_materialized_or_set_to_zero']
            for box,newbox in zip(row['true_original_phase_directed_boxes'],again['true_original_phase_directed_boxes'],strict=True):
                for endpoint in ('lower','upper'):
                    assert box[endpoint]['exact_mpf_tuple']==newbox[endpoint]['exact_mpf_tuple']
            phasechecks+=1
    frame=owner.source.frame(64,0,branch='regular')
    for phase in ((0,0),(1,1)):
        bounds=owner.primitive_bounds(frame,phase=phase)
        assert all(v.zero for orders in bounds['caps'].values() for v in orders.values())
    for call in (lambda:owner.at_candidate(dict(live),N=160),lambda:owner.at_candidate(live,N=True),
            lambda:owner.at_candidate(live,N=159),lambda:owner.at_candidate(live,N=160.0),
            lambda:owner.at_candidate(live),lambda:owner.hydrate(replace(frame)),
            lambda:owner.primitive_bounds(frame,phase=(-.1,.2)),lambda:owner.primitive_bounds(frame,phase=(.2,1.1))):
        try:call()
        except (TypeError,ValueError):rejected+=1
    assert rejected==8
    return dict(passed=True,explicit_Nminus2_and_native_unit_queries=queries,
        actual_original_common_N_radius_phase_endpoints=phasechecks,
        exact_periodic_primitive_all_slow_endpoint_traces_zero=True,runtime_candidate_source_and_phase_guards=rejected,
        original_positive_hb_sc_offset_retained=True,no_old_N7_archive_or_selected_approximate_phase=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_phase_primitive_point_evaluator_installed','sharp_O2_phase_averaging_installed',
        'actual_all_route_incoming_histories_installed','functional_terminal_identity_solved','current_whole_N_selected',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',*current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2AllNSpatialEnvelope();checks={}
    with mp.workdps(owner.ctx.dps+40):
        checks['independent_full_unpaired_exponential_remainder']=independent_full_unpaired_remainder(owner)
        print('independent_full_unpaired_exponential_remainder PASS',flush=True)
        checks['independent_actual_spatial_phase_IBP_and_full_density']=independent_nonzero_endpoint_spatial_IBP(owner)
        print('independent_actual_spatial_phase_IBP_and_full_density PASS',flush=True)
        live,checks['actual_original_source_spatial_coefficient_and_mass_binding']=actual_source_and_mass_reconstruction(owner,manifest)
        print('actual_original_source_spatial_coefficient_and_mass_binding PASS',flush=True)
        checks['actual_common_N_phase_units_and_runtime_guards']=candidate_phase_units_and_guards(owner,live,manifest)
        print('actual_common_N_phase_units_and_runtime_guards PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Complete actual original O2 all-N>=160 full changed spatial N^-2 envelopes, actual signed zero-mean primitive bounds, unpaired nonlinear bias, source/own-rate/all-endpoint/candidate phase ownership. Conservative O2 only; no sharp/all-route/terminal/global N/recursion/corrected NS claim.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Original O2 complete spatial all-N envelope PASS',flush=True);return report


if __name__=='__main__':run()
