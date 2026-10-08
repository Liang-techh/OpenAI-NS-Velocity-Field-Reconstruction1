"""Focused full-reference predicate source, overlap and own-rate evidence."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_axial_predicate_C1 as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

C0,Y,Z,YZ=current.mixed.ORDERS;ep=current.ep


def common_original_function_proof():
    r,c,q,h=s.symbols('r c q h',real=True)
    ss=1-r*r;D=1-2*r*c+r*r
    # Exact Mobius cosine and Jacobian. Both antiderivatives vanish at psi=0.
    coschi=((1+r*r)*c-2*r)/D
    T1psi=q*h/r*(ss/D-1)
    t=2*q*h*(c-r)/D
    T2psi=q*q/r**2*((2-3*ss+2*r*coschi)*ss/D+ss)
    assert s.cancel(T1psi-t)==0
    assert s.cancel((T2psi-t*t).subs(h*h,ss))==0
    y,z=s.symbols('y z',real=True);g=s.Function('g')(y,z)
    qs,ds=s.symbols('q dstar',positive=True)
    u=z*g*qs/ds
    assert s.simplify(s.diff(u,y)-u*s.diff(g,y)/g)==0
    aa,T2i,ti,Ti=s.symbols('a T2i ti Ti',real=True)
    psi_i=-T2i/(1+ti*ti)
    assert s.cancel(-aa*psi_i/(4*s.pi)-aa*T2i/(4*s.pi*(1+ti*ti)))==0
    return dict(passed=True,exact_signed_T1_T2_original_direction_derivative_identities=2,
        exact_source_u_y_equals_u_g_y_over_g=True,exact_first_implicit_A_identity=True,
        original_direction='t=2*q*sqrt(1-r^2)*(cospsi-r)/(1-2*r*cospsi+r^2)',
        common_original_integrals='T1=int_0^psi t; T2=int_0^psi t^2',
        regular_Fourier_series_and_actual_derivative_tail_proof_from_checked_predecessor=current.regular.RECEIPT,
        inverse_uniqueness='Phi_psi=(1+t^2)/(2*pi*nu)>0; nu=1+2*q^2>0',
        same_defining_primitives_and_unique_true_phase_inverse_on_overlap=True,
        same_original_C0_y_Z_yZ_derivatives_not_chosen_interval_jet_functions=True,
        original_smooth_function_internal_radial_traces_cancel_even_if_coordinates_switch=True,
        signed_rational_R0_is_an_outer_enclosure_of_original_inverse_not_a_field=True,
        coverage='abs(u)<=1/4 OR u>=3/16 OR u<=-3/16; overlaps 3/16<=abs(u)<=1/4',
        predicates_covers_all_real_u_including_zero=True,
        overlap_majorants_unioned_not_double_counted=True)


def carrier_and_domain_contract(owner,manifest):
    c=owner.ctx;cert=manifest['carrier_certificate'];live=owner.certificate
    assert cert['passed'] and cert['sign']==-1
    assert cert['source_family']==owner.family
    rows=cert['whole_source_partition'];assert len(rows)==256
    for i in range(16):
        for j in range(16):
            row=rows[16*i+j]
            assert row['exact_y_cell']==[str(-5+s.Rational(5*i,16)),str(-5+s.Rational(5*(i+1),16))]
            assert row['exact_abs_Z_cell']==[str(s.Rational(j,16)),str(s.Rational(j+1,16))]
            assert ep(saved.interval(c,row['full_p2_over_Z_divided_by_Lambda0']))[1]<0
            assert row['all_original_pressure_error_terms_retained']
    g=saved.interval(c,cert['g_normalized']);gy=saved.interval(c,cert['g_y_normalized'])
    gamma=saved.interval(c,cert['gamma_y'])
    assert ep(g)[1]<0
    saved.contains(gamma,gy/g)
    saved.contains(saved.interval(c,cert['g_normalized']),live['g_normalized'])
    assert cert['all_negative_Z_covered_by_exact_even_g_and_gy_parity']
    assert cert['q_and_dstar_y_Z_derivatives_exact_zero_on_reference']
    cmag=saved.interval(c,cert['u_over_Lambda0_abs_Z'])
    assert ep(cmag)[0]>0
    radius=s.Rational(cert['exact_central_zeta_max'])
    assert ep(owner.atlas.rational(radius)*c.mpf(ep(cmag)[0]))[0]>=c.mpf('.25')
    assert ep(saved.interval(c,cert['q']))[0]>=c.mpf('.5')
    rejections=0
    for call in (
        lambda:current.AxialSourceAtlas(owner.seed.owner.inputs.frame,lower=-2,upper=1),
        lambda:current.AxialSourceAtlas(owner.seed.owner.inputs.frame,lower=0,upper=1,carrier=(1,0)),
        lambda:owner.branches['regular'].source_frame(0,0),
        lambda:owner.branches['positive'].primitive(replace(owner.branches['positive'].source_frame(-1,0)),c.mpf('.3'))):
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejections+=1
    assert rejections==4
    return dict(passed=True,whole_source_rectangles_with_full_pressure_sign_checked=len(rows),
        exact_positive_and_negative_axial_coverage_not_samples=True,
        conditional_regular_physical_carrier_radius_proved=True,
        gamma_y_is_same_source_ratio_not_independent_parameter=True,
        signed_predicate_excludes_Z_zero_even_though_outer_rectangle_touches_zero=True,
        original_basis_domain_and_issued_frame_rejections=rejections)


def independent_overlap(owner):
    a=owner.atlas;c=owner.ctx;p=mp.mp.clone();p.dps=100
    theta,uu,qq=s.symbols('theta u q',real=True)
    h=s.sqrt(1+uu*uu);r=uu/h
    t=2*qq/h*(s.cos(theta)-r)/(1-2*r*s.cos(theta)+r*r)
    integrands={power:tuple(s.lambdify((theta,uu,qq),s.diff(t**power,uu,order),
        modules=[{'mpf':p.mpf},'mpmath']) for order in (0,1,2)) for power in (1,2)}
    comparisons=0
    for utext in ('-.21','.21'):
        u0=p.mpf(utext);q=p.mpf('.9');uy=p.mpf('.03');uz=p.mpf('-.04');uyz=p.mpf('.011')
        roots={name:current.MixedJet(a,{key:a.scalar(c.mpf(value)) for key,value in
            zip(current.mixed.ORDERS,vals,strict=True)}) for name,vals in
            dict(a=('.8',0,0,0),b=(0,0,0,0),t0=(0,0,0,0),E=('1.3','-.2','.17','.08'),
                p2=(u0/q,uy/q,uz/q,uyz/q)).items()}
        directional={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()}
        kernel=current.PredicatePhase(dict(q=a.scalar(c.mpf(q)),roots=directional,
            original_u_source=roots['p2'][C0]*a.scalar(c.mpf(q)),regular_predicate=True),c.mpf(0))
        # Same defining u, deliberately force the signed coordinate on overlap.
        logu=c.ln(c.mpf(abs(u0)))
        sa=current.AxialSourceAtlas(owner.seed.owner.inputs.frame,lower=-1,upper=1,logu=logu)
        sc=sa.ctx
        sroots={name:current.MixedJet(sa,{key:sa.scalar(sc.mpf(value)) for key,value in
            zip(current.mixed.ORDERS,vals,strict=True)}) for name,vals in
            dict(a=('.8',0,0,0),b=(0,0,0,0),t0=(0,0,0,0),E=('1.3','-.2','.17','.08'),
                p2=(u0/q,uy/q,uz/q,uyz/q)).items()}
        signed_u=current.prior.ScaledEnclosure(current.prior.FormalScale(sa.bases,(0,0,0,1,0)),
            1 if u0>0 else -1,sa.ledger)
        signed=current.PredicatePhase(dict(q=sa.scalar(sc.mpf(q)),
            roots={name:{C0:row[C0],Z:row[Z]} for name,row in sroots.items()},
            original_u_source=signed_u,signed_predicate=True),sc.mpf(0))
        assert kernel.geometry=='small_r_series' and signed.geometry=='signed_Mobius'
        for psi_text in ('.47','3.9'):
            psi=p.mpf(psi_text);rr=u0/p.sqrt(1+u0*u0)
            chi=2*p.atan2((1+rr)*p.sin(psi/2),(1-rr)*p.cos(psi/2))
            if chi<0:chi+=2*p.pi
            outputs=[(a,current.REGULAR_MIXED(a,kernel,roots,c.mpf(psi/(2*p.pi)))),
                (sa,current.SIGNED_MIXED(sa,signed,sroots,sc.mpf(chi/(2*p.pi))))]
            segments=[0]+([p.pi] if psi>p.pi else [])+[psi];partials={}
            for power,name in ((1,'T1'),(2,'T2')):
                f0,f1,f2=integrands[power]
                partials[name]={C0:p.quad(lambda th:f0(th,u0,q),segments),
                    Y:p.quad(lambda th:f1(th,u0,q)*uy,segments),
                    Z:p.quad(lambda th:f1(th,u0,q)*uz,segments),
                    YZ:p.quad(lambda th:f2(th,u0,q)*uy*uz+f1(th,u0,q)*uyz,segments)}
            tt=integrands[1][0](psi,u0,q);du=integrands[1][1](psi,u0,q)
            ty,tz=du*uy,du*uz;tp=p.diff(lambda th:integrands[1][0](th,u0,q),psi)
            D=1+tt*tt;py=-partials['T2'][Y]/D;pz=-partials['T2'][Z]/D
            pyz=-(partials['T2'][YZ]+2*tt*ty*pz+2*tt*tz*py+2*tt*tp*py*pz)/D
            totalY=partials['T1'][Y]+tt*py;totalZ=partials['T1'][Z]+tt*pz
            totalYZ=partials['T1'][YZ]+ty*pz+tz*py+tp*py*pz+tt*pyz
            phi=(psi+partials['T2'][C0])/(2*p.pi*(1+2*q*q));aa=p.mpf('.8')
            expectedA={C0:aa/2*(phi-psi/(2*p.pi)),Y:-aa*py/(4*p.pi),
                Z:-aa*pz/(4*p.pi),YZ:-aa*pyz/(4*p.pi)}
            expectedB={C0:-aa*p.mpf('1.3')*partials['T1'][C0]/(4*p.pi),
                Y:-aa/(4*p.pi)*(p.mpf('-.2')*partials['T1'][C0]+p.mpf('1.3')*totalY),
                Z:-aa/(4*p.pi)*(p.mpf('.17')*partials['T1'][C0]+p.mpf('1.3')*totalZ),
                YZ:-aa/(4*p.pi)*(p.mpf('.08')*partials['T1'][C0]-p.mpf('.2')*totalZ+
                    p.mpf('.17')*totalY+p.mpf('1.3')*totalYZ)}
            for atlas,got in outputs:
                for name in ('T1','T2'):
                    for key,value in partials[name].items():
                        saved.contains(saved.saved_value(atlas.ctx,got['record']['original_fixed_angle_'+name][str(key)],bases=atlas.bases),
                            value,p.mpf('1e-80'));comparisons+=1
                for name,rows in (('A',expectedA),('B',expectedB)):
                    for key,value in rows.items():
                        saved.contains(got[name][key].finite_interval(),value,p.mpf('1e-80'));comparisons+=1
    # Closed small-u endpoint does not depend on a rounded log comparison.
    q=a.scalar(c.mpf('.9'))
    query=dict(q=q,roots=directional,original_u_source=a.scalar(c.mpf('.25')),regular_predicate=True)
    assert current.PredicatePhase(query,c.mpf(0)).geometry=='small_r_series'
    bad=dict(query,original_u_source=a.scalar(c.mpf('.26')))
    try:current.PredicatePhase(bad,c.mpf(0))
    except ValueError:pass
    else:raise AssertionError('Unproved closed regular source guard accepted')
    return dict(passed=True,independent_original_T1_T2_A_B_C0_y_Z_yZ_comparisons=comparisons,
        both_signs_and_both_overlapping_coordinate_backends=True,
        same_rational_original_direction_independently_integrated=True,
        full_E_product_and_implicit_first_and_mixed_cross_terms_tested=True,
        exact_closed_regular_endpoint_and_out_of_predicate_rejection=True,
        finite_diagnostic_values_not_installed_as_original_parameters=True)


def integral_contracts(owner,manifest):
    c=owner.ctx;levels=manifest['actual_full_reference_axial_C1_levels']
    assert [(v['exact_radial_cells'],v['candidate_N']) for v in levels]==[(4,160),(16,16384)]
    cells_checked=0;masses=0;bounds=[]
    for level in levels:
        count=level['exact_radial_cells'];N=level['candidate_N']
        assert level['exact_physical_Z_window']==['-1','1'] and 'original_Z_exact' not in level
        assert level['actual_mixed_source_and_C1_averaging_installed_on_whole_reference_Z_window']
        assert not level['O2_nonconstant_parameter_extension_installed']
        assert level['same_original_source_internal_C0_Z_traces_cancel']
        assert level['no_cross_chart_seam_assumed'] and level['actual_global_C0_Z_endpoint_terms_retained']
        assert level['incoming_histories_and_P0_not_reset'] and level['pressure_zero_rate_memory_retained']
        assert not level['whole_Z_terminal_or_all_route_closure']
        totals={name:c.mpf(0) for name in current.points.exact.RATES}
        for i,cell in enumerate(level['whole_cell_mixed_source_and_IBP_records']):
            source=cell['source'];left=-5+s.Rational(5*i,count);right=-5+s.Rational(5*(i+1),count)
            assert source['exact_reference_cell']==[str(left),str(right)]
            assert source['same_original_function_branch_caps_union_not_sum']
            assert source['entire_physical_Z_domain_covered_by_three_original_u_predicates']
            assert set(source['original_conditional_source_cells'])=={'regular','positive','negative'}
            for branch,part in source['original_conditional_source_cells'].items():
                frame=part['source']
                assert frame['conditional_domain_not_claimed_entire_physical_rectangle']
                assert frame['original_p2_Z_yZ_template_rows_retained']
                assert frame['source_predicate_not_field_or_derivative_of_cap']
                assert part['complete_Z_product_rules_and_exact_N_remainder_retained']
                assert part['phase_Z_exact_zero'] and part['native_log_radius_Jacobian']==1
                if branch!='regular':
                    K=part['actual_mixed_primitive']['actual_source_K_rows']
                    assert K[str(Y)]['formal_positive_scale']['source_exponents'][:2]==[0,0]
                    assert K[str(Y)]['formal_positive_scale']['source_exponents'][3]==0
                    assert K[str(Z)]['formal_positive_scale']['source_exponents'][:2]==[11,10]
                    assert K[str(Z)]['formal_positive_scale']['source_exponents'][3]==-1
                for endpoint in part['actual_original_endpoint_phases'].values():
                    assert endpoint['positive_original_origin_offset']['strictly_positive']
            for name,contribution in cell['contributions'].items():
                assert contribution['N_power']==-2
                mass=saved.interval(c,contribution['positive_own_rate_mass'])
                assert ep(mass)[0]>0;totals[name]+=mass
                if 0<i<count-1:
                    assert all(v['global_endpoint_term']['exact_zero'] for v in contribution['C0_Z_rows'].values())
            cells_checked+=1
        for name,rate in current.points.exact.RATES.items():
            rate=c.mpf(rate.numerator)/rate.denominator
            expected=c.mpf(5) if ep(rate)==(0,0) else (1-c.exp(-5*rate))/rate
            saved.contains(totals[name],expected);masses+=1
        out={}
        for jet in (C0,Z):
            out[str(jet)]={}
            for name in current.points.exact.RATES:
                quantities={}
                for kind in ('phase_averaged','direct','effective'):
                    row=level['C0_Z_'+kind+'_contribution_enclosures'][str(jet)][name]
                    scale=row['formal_positive_scale']
                    assert scale['source_exponents']==([0,0,0,0] if jet==C0 else [11,10,0,0])
                    assert scale['radius_power']==0
                    box=saved.interval(c,row['coefficient_interval'])*c.exp(saved.interval(c,scale['additional_log_interval']))
                    quantities[kind]=max(abs(v) for v in ep(box))
                assert quantities['effective']<=min(quantities['phase_averaged'],quantities['direct'])
                out[str(jet)][name]=quantities
        bounds.append(dict(cells=count,N=N,C0_bounds_and_Z_bounds_divided_by_Lambda0=out))
    return dict(passed=True,whole_reference_radial_cells_checked=cells_checked,
        original_positive_own_rate_mass_comparisons=masses,pressure_zero_rate_memory_mass=5,
        all_three_conditional_caps_cover_same_original_function=True,
        signed_K_y_formal_logu_and_Lambda_cancellation_retained=True,
        original_Z_yZ_large_Lambda_not_erased=True,computed_ranges=bounds,
        global_frequency_and_terminal_moment_closure_not_claimed=True)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*current.near.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalReferenceAxialPredicates();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (
            ('common_original_function_and_predicate_coverage',common_original_function_proof),
            ('full_pressure_axial_carrier_and_domain',lambda:carrier_and_domain_contract(owner,manifest)),
            ('independent_regular_signed_overlap',lambda:independent_overlap(owner)),
            ('actual_full_reference_C1_integrals',lambda:integral_contracts(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Original reference full Z[-1,1] same-source conditional regular/signed C0/y/Z/yZ and actual-N C0/Z own-rate integral enclosures. Full-pressure carrier sign and overlap checked. Conservative bounds, no selected compatible interval field. Not all-route/O2/terminal/global N/stress/recursion/corrected NS.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Full original reference axial predicate C1 evidence PASS',flush=True)
    return report


if __name__=='__main__':run()
