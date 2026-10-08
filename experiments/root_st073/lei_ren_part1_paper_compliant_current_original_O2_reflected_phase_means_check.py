"""Same-function pair, true phase measure and actual mean contribution check."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_reflected_phase_means as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;C0,Y,Z,YZ=current.mixed.ORDERS


def independent_pair_mixed_fixtures(owner):
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas;c=owner.ctx
    p=mp.mp.clone();p.dps=170;comparisons=0;cases=0
    for N in (160,257):
        for Vtext in ('-.21','.21'):
            for Atext in ('-.7','0','1e-50','.6'):
                rows=dict(E=('1.1','.2','-.3','.04'),V=(Vtext,'.07','.4','-.05'),
                    A=(Atext,'.17','-.13','.11'),B=('.53','-.09','.23','-.08'))
                jets={name:current.MixedJet(a,{order:a.scalar(c.mpf(value)) for order,value in zip(current.mixed.ORDERS,values,strict=True)})
                    for name,values in rows.items()}
                actual=current.reflected_pair_densities(a,**jets,N=N);yy=p.mpf('.41');zz=p.mpf('-.37')
                def function(name,y,z):
                    v,dy,dz,dyz=map(p.mpf,rows[name]);return v+dy*(y-yy)+dz*(z-zz)+dyz*(y-yy)*(z-zz)
                def densities(E,V):return dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
                def pair(key,y,z):
                    E,V,A,B=(function(name,y,z) for name in ('E','V','A','B'))
                    old=densities(E,V)
                    plus=densities(E*p.exp(A/N),V+B/N);minus=densities(E*p.exp(-A/N),V-B/N)
                    return (plus[key]+minus[key])/2-old[key]
                for key,jet in actual.items():
                    targets={C0:pair(key,yy,zz),Y:p.diff(lambda y:pair(key,y,zz),yy),Z:p.diff(lambda z:pair(key,yy,z),zz),
                        YZ:p.diff(lambda y:p.diff(lambda z:pair(key,y,z),zz),yy)}
                    for order,target in targets.items():
                        saved.contains(jet[order].finite_interval(),target,p.mpf('1e-120'));comparisons+=1
                assert all(actual['m'][order].zero for order in current.mixed.ORDERS);cases+=1
    return dict(passed=True,independent_original_exponential_pair_C0_y_Z_yZ_comparisons=comparisons,
        compatible_actual_mixed_function_fixtures=cases,nonzero_V_term_and_all_cross_terms_included=True,
        tiny_and_zero_A_stable_hyperbolic_chain_included=True)


def independent_original_true_phase_measure():
    """Finite diagnostics integrate the original angle loop and true measure.

    These finite parameters are independent diagnostic cases; they do not
    replace the original extraordinarily scaled source or prove admission.
    """
    p=mp.mp.clone();p.dps=65;cases=0;comparisons=0;measurechecks=0
    for utext in ('-.22','.22'):
        for qtext in ('.2','.001'):
            u,q=p.mpf(utext),p.mpf(qtext);r=u/p.sqrt(1+u*u);h=1/p.sqrt(1+u*u);ss=h*h
            a=p.mpf('1.2');E=p.mpf('.8');V=p.mpf('.21') if u>0 else p.mpf('-.21');N=160;nu=1+2*q*q
            def roots(psi):
                chi=2*p.atan2((1+r)/(1-r)*p.sin(psi/2),p.cos(psi/2))
                T1=q*h/r*(chi-psi)
                T2=q*q/(r*r)*((2-3*ss)*chi+ss*psi+2*r*p.sin(chi))
                phi=(psi+T2)/(2*p.pi*nu);A=a*(phi-psi/(2*p.pi))/2;B=-a*E*T1/(4*p.pi)
                t=2*q*h*(p.cos(psi)-r)/(1-2*r*p.cos(psi)+r*r)
                return A,B,(1+t*t)/(2*p.pi*nu),t
            full=[0,p.pi/2,p.pi,3*p.pi/2,2*p.pi];half=[0,p.pi/2,p.pi]
            weight=p.quad(lambda psi:roots(psi)[2],full)
            tmean=p.quad(lambda psi:roots(psi)[3],full)
            t2mean=p.quad(lambda psi:roots(psi)[3]**2,full)
            assert abs(weight-1)<p.mpf('1e-55') and abs(tmean)<p.mpf('1e-55')
            assert abs(t2mean-4*p.pi*q*q)<p.mpf('1e-55');measurechecks+=3
            def direct(key,psi):
                A,B,w,t=roots(psi);de=E*p.expm1(A/N);dv=B/N
                return dict(m=dv,h=de,k=V*de+E*dv+de*dv,
                    e=2*V*dv+dv*dv-E*de-de*de/2,p=E*de+de*de/2)[key]*w
            def paired(key,psi):
                A,B,w,t=roots(psi);x,z=A/N,B/N
                return dict(m=p.mpf(0),h=E*(p.cosh(x)-1),k=V*E*(p.cosh(x)-1)+E*z*p.sinh(x),
                    e=z*z-E*E*(p.cosh(2*x)-1)/2,p=E*E*(p.cosh(2*x)-1)/2)[key]*w*2
            means={}
            for key in current.integrals.KEYS:
                original=p.quad(lambda psi:direct(key,psi),full)
                pair=p.quad(lambda psi:paired(key,psi),half)
                assert abs(original-pair)<p.mpf('1e-55'),(utext,qtext,key);comparisons+=1;means[key]=pair
            assert means['m']==0 and means['h']>0 and means['p']>0 and means['e']+means['p']>0
            cases+=1
    return dict(passed=True,independent_original_loop_phase_mean_comparisons=comparisons,
        independent_phase_measure_and_original_period_moment_checks=measurechecks,finite_positive_q_both_sign_cases=cases,
        integrated_true_phase_measure_not_uniform_free_angle=True,
        pair_integrands_integrated_before_claiming_full_means=True,
        finite_parameters_not_original_source_replacements=True)


def actual_mean_domain_and_mass(owner,manifest):
    c=owner.ctx;p=mp.mp.clone();p.dps=150
    record=manifest['actual_original_O2_true_phase_means_and_finite_N_bias'];rows=record['original_source_true_phase_mean_records']
    assert record['source_family']==owner.family and record['explicit_candidate_N']==160
    assert record['exact_y_window']==['0','1'] and record['exact_Z_range']==['-1','1']
    assert record['actual_true_phase_interval_inverse_queries']==768 and len(rows)==64
    assert record['own_rates']==current.integrals.five.RATES and record['normalized_own_units']==current.integrals.five.UNITS
    assert record['mean_Duhamel_contribution_not_complete_oscillatory_contribution']
    assert record['oscillatory_spatial_density_integral_remainder_enclosed'] is False
    totals=[{key:c.mpf(0) for key in current.integrals.KEYS} for unused in range(2)]
    pairrows=0;means=0;exports=0;masschecks=0;inversechecks=0
    for index,row in enumerate(rows):
        assert row['exact_y_cell']==[str(s.Rational(index,64)),str(s.Rational(index+1,64))]
        assert row['branch_overlap_hulled_not_summed'] and row['phase_reflection_weight_applied_once']
        assert row['no_extra_R_Jacobian_or_period_count']
        for branchrecord in row['conditional_source_true_phase_mean_records']:
            branch=branchrecord['branch'];frame=owner.source.frame(64,index,branch=branch);a=frame.roots['q'].atlas
            assert branchrecord['actual_original_predicate']==frame.record['actual_original_source_predicate']
            assert branchrecord['exact_outer_Z_bounds']==frame.record['exact_outer_Z_bounds']
            parts=branchrecord['complete_original_half_phase_partition'];assert len(parts)==4
            mean={key:current.MixedJet.constant(a,0) for key in current.integrals.KEYS}
            for j,part in enumerate(parts):
                left,right=s.Rational(j,8),s.Rational(j+1,8)
                assert part['exact_true_phase_cell']==[str(left),str(right)]
                assert part['exact_reflected_true_phase_cell']==[str(1-right),str(1-left)]
                assert part['exact_pair_to_full_mean_weight']=='1/4'
                assert part['source_family']==owner.family and part['source_level']==64 and part['source_index']==index
                assert part['paired_same_original_function_not_two_independent_caps']
                assert part['nonzero_original_V_term_in_k_retained']
                inverse=part['selected_original_interval_inverse']
                assert inverse['chart'] in ('psi','E') and 'directed endpoint' in inverse['bracket_proof']
                saved.contains(saved.interval(c,inverse['phase_image']),saved.interval(c,part['original_phase_interval']),mp.mpf('1e-200'));inversechecks+=1
                for key,jetrows in part['actual_reflected_pair_density_C0_y_Z_yZ'].items():
                    jet=current.MixedJet(a,{order:current.integrals.restore_value(a,jetrows[str(order)]) for order in current.mixed.ORDERS})
                    for order in current.mixed.ORDERS:
                        assert jetrows[str(order)]['encloses_original_source_function'] and not jetrows[str(order)]['point_value_selected'];pairrows+=1
                    if key=='m':assert all(jet[order].zero for order in current.mixed.ORDERS)
                    mean[key]=mean[key]+jet*(c.mpf(1)/4)
            for key,jet in mean.items():
                for order in current.mixed.ORDERS:
                    stored=current.integrals.restore_value(a,branchrecord['actual_true_phase_mean_C0_y_Z_yZ'][key][str(order)])
                    assert stored.scale.powers==jet[order].scale.powers
                    assert ep(stored.scale.offset)==ep(jet[order].scale.offset) and ep(stored.coefficient)==ep(jet[order].coefficient);means+=1
                for j,order in enumerate((C0,Z)):
                    value=current.integrals.export_range(a,jet[order],derivative_unit=order==Z)
                    exported=saved.interval(c,branchrecord['normalized_mean_C0_Z_over_source_unit'][j][key])
                    saved.contains(exported,value,mp.mpf('1e-200'))
                    saved.contains(saved.interval(c,row['normalized_mean_C0_Z_unions'][j][key]),value,mp.mpf('1e-200'));exports+=1
                if key in ('h','p'):assert ep(current.integrals.export_range(a,jet[C0],derivative_unit=False))[0]>=0
        left,right=p.mpf(index)/64,p.mpf(index+1)/64
        for key,rate in current.integrals.five.RATES.items():
            rr=p.mpf(rate);actual=right-left if not rr else (p.exp(-rr*(1-right))-p.exp(-rr*(1-left)))/rr
            mass=saved.interval(c,row['positive_own_rate_final_endpoint_masses'][key]);assert ep(mass)[0]>0
            saved.contains(mass,actual,p.mpf('1e-120'));masschecks+=1
            for j in range(2):
                contribution=saved.interval(c,row['normalized_mean_C0_Z_unions'][j][key])*mass
                saved.contains(saved.interval(c,row['normalized_mean_Duhamel_contributions'][j][key]),contribution,mp.mpf('1e-200'))
                totals[j][key]+=contribution
    for j in range(2):
        for key in current.integrals.KEYS:saved.contains(saved.interval(c,record['all_normalized_C0_Z_mean_bias_contributions'][j][key]),totals[j][key],mp.mpf('1e-200'))
    assert record['true_phase_mean_Duhamel_C0']['m']['exact_zero'] and record['true_phase_mean_Duhamel_ordinary_Z']['m']['exact_zero']
    return dict(passed=True,original_true_phase_interval_inverses=inversechecks,native_original_reflected_pair_density_rows=pairrows,
        exact_native_four_slow_mean_rows_reconstructed=means,original_C0_Z_normalized_exports_and_branch_unions=exports,
        independent_positive_own_rate_masses=masschecks,continuous_original_y_Z_predicate_cover=True,
        nonlinear_mean_bias_not_complete_spatial_oscillatory_contribution=True)


def source_and_hyperbolic_guards(owner):
    from dataclasses import replace
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas;rejected=0
    for call in (lambda:owner.pair(replace(frame),0,s.Rational(1,8),N=160),
            lambda:owner.pair(frame,0,1,N=160),lambda:owner.pair(frame,0,0,N=160),
            lambda:owner.pair(frame,0,s.Rational(1,8),N=True),
            lambda:current.stable_hyperbolic_changes(a,current.MixedJet.constant(a,a.scalar(2)))):
        try:call()
        except ValueError:rejected+=1
    assert rejected==5
    return dict(passed=True,copied_owner_true_half_phase_candidate_N_and_unbounded_hyperbolic_guards=rejected)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('sharp_O2_phase_averaging_installed','oscillatory_spatial_integral_remainder_enclosed','actual_all_route_incoming_histories_installed',
        'functional_terminal_identity_solved','current_whole_N_selected','all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2ReflectedMeans();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('same_original_phase_reflection_and_pair_calculus',current.reflection_theorem),
                ('independent_nonlinear_pair_mixed_functions',lambda:independent_pair_mixed_fixtures(owner)),
                ('independent_original_loop_true_phase_measure',independent_original_true_phase_measure),
                ('actual_continuous_source_phase_means_and_Duhamel_bias',lambda:actual_mean_domain_and_mass(owner,manifest)),
                ('source_phase_and_hyperbolic_guards',lambda:source_and_hyperbolic_guards(owner))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Actual source true-phase reflection, stable mixed nonlinear pair, full original phase mean C0/Z and mean Duhamel contribution. Finite-N biases retained, spatial oscillatory remainder and sharp averaging/terminal/global N/all-route/full reconstruction still open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Original O2 true phase means and finite-N bias PASS',flush=True);return report


if __name__=='__main__':run()
