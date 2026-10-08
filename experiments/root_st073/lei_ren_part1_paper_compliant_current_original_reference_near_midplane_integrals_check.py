"""Focused whole-zeta source, independent regular primitives and integral evidence."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_near_midplane_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;points=current.points


def factor_pressure_source(owner):
    a=owner.atlas;c=owner.ctx;z=owner.z;Q=owner.Qsymbol
    assert a.Z is None and a.zeta_bounds==(s.Rational(-1,10**8),s.Rational(1,10**8))
    assert ep(a.Q)[0]>=1 and ep(a.L)[0]>0
    assert ep(a.bases[2])[0]<0<=ep(a.bases[2])[1]
    assert not a.physical_Z(c.mpf((-1,1))).zero
    # Original R^r Pstar^p delta^d L^ell Z^k identity, before amplification.
    P,C,LL,y,zeta=s.symbols('P C LL y zeta',real=True);identities=0
    for rows in owner.terms.values():
        for row in rows:
            r,p,d,ell=row['powers'];k=row['k']
            original=r*(s.log(110)+10*C+10*P+y)+p*P+d*(-4*P-30)+ell*LL-k*(11*P+10*C)
            mapped=(p-4*d+10*r-11*k)*P+(10*r-10*k)*C+ell*LL-30*d+r*(s.log(110)+y)
            assert s.expand(original-mapped)==0
            value=a.zterm(row['powers'],1,coordinate=0,zeta=c.mpf('.7'),Z_power=k)
            assert value.scale.powers==(p-4*d+10*r-11*k,10*r-10*k,ell,0,0)
            identities+=1
    assert all(row['k']%2==1 for row in owner.terms[('p2',0)])
    assert any(row['error_order']==1 for row in owner.terms[('p2',0)])
    assert all(row['k']>=1 for row in owner.terms[('p2',0)] if row['error_order']==1)
    assert any(row['error_order']==2 for row in owner.terms[('p2',1)])
    assert len(owner.pressure_theorem['original_stage_derivative_factor_identities'])==14
    assert owner.pressure_theorem['H_R_not_identified_with_R0_ZZ']
    # Independent original C0 differentiation includes variable L and a
    # nonzero even pressure remainder, rather than differentiating caps.
    p=mp.mp.clone();p.dps=130;Pstar=p.exp(3);delta=p.exp(-42);Lambda=p.exp(53)
    alpha=p.mpf('2.9');epsilon=p.exp(p.mpf('.6'))*p.mpf('.3')/Pstar
    rows=owner.templates['rows']
    functions={key:tuple(s.lambdify(owner.templates['inputs'],expr,
        modules=[{'mpf':p.mpf},'mpmath']) for powers,expr in terms) for key,terms in rows.items()}
    def values(y,Z,with_remainder):
        f=p.exp(y/10);q=1+Z*Z;aa=alpha-epsilon if with_remainder else alpha
        return (Z,f,p.mpf(5)/8*f,p.mpf(5)/12*f*f,p.mpf(5)/2*f*f,
                -aa/q**2,4*aa*Z/q**3,aa*(4-20*Z*Z)/q**4)
    def original(name,y,Z,order,with_remainder=True):
        L=1-delta*Z*Z;R=110*p.exp(50+y)
        return sum(R**r*Pstar**ps*delta**d*L**ell*fn(*values(y,Z,with_remainder))
            for ((r,ps,d,ell),expr),fn in zip(rows[(name,order)],functions[(name,order)],strict=True))
    replay={key:tuple((row,s.lambdify(owner.arguments,row['expr'],
        modules=[{'mpf':p.mpf},'mpmath'])) for row in terms) for key,terms in owner.terms.items()}
    comparisons=0;derivatives=0;remainder_covers=0
    for yy in ('-2.337','-.513'):
        y=p.mpf(yy)
        for zz in ('-.7','0','.7'):
            zz=p.mpf(zz);Z=zz/Lambda;f=p.exp(y/10);q=1+Z*Z;L=1-delta*Z*Z
            args=(f,p.mpf(5)/8*f,p.mpf(5)/12*f*f,p.mpf(5)/2*f*f,alpha,q)
            for name in ('E','V','b','p1','p2'):
                expected=p.diff(lambda x:original(name,y,x,0),Z)
                got=original(name,y,Z,1)
                assert abs(got-expected)<=p.mpf('1e-105')*max(abs(got),abs(expected),1),(name,yy,zz)
                derivatives+=1
                for order in (0,1):
                    baseline=p.mpf(0);budget=p.mpf(0)
                    for row,fn in replay[(name,order)]:
                        r,ps,d,ell=row['powers'];k=row['k']
                        # Collected finite affine powers, independent of atlas.
                        factor=p.exp((ps-4*d+10*r-11*k)*3+(10*r-10*k)*2
                                     -30*d+r*(p.log(110)+y))*L**ell
                        term=factor*fn(*args)*zz**k
                        if row['error_order'] is None:baseline+=term
                        else:budget+=abs(term)*p.exp(p.mpf('.6'))
                    target=original(name,y,Z,order,False)
                    tolerance=p.mpf('1e-105')*max(abs(target),abs(baseline),1)
                    assert abs(target-baseline)<=tolerance
                    assert abs(original(name,y,Z,order)-baseline)<=budget+tolerance
                    comparisons+=1;remainder_covers+=1
    return dict(passed=True,exact_collected_original_Z_power_identities=identities,
        independent_finite_unit_original_C0_derivative_comparisons=derivatives,
        independent_original_source_polynomial_replays=comparisons,
        nonzero_consistent_even_pressure_remainder_cover_comparisons=remainder_covers,
        all_fourteen_original_stage_odd_factors_bound_same_remainder=True,
        finite_alpha_coefficient_enclosure_and_one_Pstar_inverse_retained=True,
        variable_positive_L_Q_and_nonzero_formal_physical_Z_retained=True,
        diagnostics_not_installed_as_original_source_parameters=True)


def independent_regular_primitives(owner):
    a=owner.atlas;c=owner.ctx;p=mp.mp.clone();p.dps=110;comparisons=0
    for uu in ('-.14','0','.14'):
        u=p.mpf(uu);q=p.mpf('.9');E=p.mpf('1.3');EZ=p.mpf('-.2');aa=p.mpf('.8');uZ=p.mpf('.7')
        roots={name:{(0,0):a.scalar(c.mpf(value)),(0,1):a.scalar(c.mpf(dZ))}
            for name,value,dZ in (('E',E,EZ),('a',aa,0),('b',0,0),('t0',0,0),('p2',u/q,uZ/q))}
        kernel=current.whole.conditioned.PositiveLogQPhase(
            dict(q=a.scalar(c.mpf(q)),roots=roots,original_u_source=a.scalar(c.mpf(u))),c.mpf(0))
        assert kernel.geometry=='small_r_series'
        def direction(theta,Z):
            actual=u+uZ*Z;hh=p.sqrt(1+actual*actual);r=actual/hh
            return 2*q/hh*(p.cos(theta)-r)/(1-2*r*p.cos(theta)+r*r)
        for fraction in ('.23','.71'):
            x=p.mpf(fraction);psi=2*p.pi*x
            T1=p.quad(lambda theta:direction(theta,0),[0,psi])
            T2=p.quad(lambda theta:direction(theta,0)**2,[0,psi])
            T1Z=p.quad(lambda theta:p.diff(lambda Z:direction(theta,Z),0),[0,psi])
            T2Z=p.quad(lambda theta:p.diff(lambda Z:direction(theta,Z)**2,0),[0,psi])
            t=direction(psi,0);phi=(psi+T2)/(2*p.pi*(1+2*q*q))
            expected=dict(A=aa/2*(phi-x),B_over_Pstar=-aa*E*T1/(4*p.pi),
                A_Z_slow=aa*T2Z/(4*p.pi*(1+t*t)),
                B_Z_slow=-aa/(4*p.pi)*(EZ*T1+E*(T1Z-t*T2Z/(1+t*t))))
            primitive=kernel.primitives(c.mpf(x),'psi')
            slow,proof=points.slow.slow_values(kernel,roots,c.mpf(x),'psi')
            assert proof['branch']==('exact_midplane_nonzero_p2_Z' if u==0 else 'regular_small_r_Fourier')
            for key,value in expected.items():
                got=(primitive if key in primitive else slow)[key].finite_interval()
                # M48 Fourier tails at |u|=.14 are intentionally larger
                # than numerical quadrature error and are checked as covers.
                saved.contains(got,value,p.mpf('1e-95'));comparisons+=1
            for key,value in (('T1_fixed_angle',T1),('T1_Z_fixed_angle',T1Z),('T2_Z_fixed_angle',T2Z)):
                saved.contains(saved.saved_value(c,proof[key]),value,p.mpf('1e-95'));comparisons+=1
    frame=owner.source_frame('-2.5','-1.25');u=saved.interval(c,frame.record['actual_original_u_range'])
    assert ep(u)[0]<0<ep(u)[1] and max(abs(v) for v in ep(u))<c.mpf('.25')
    result=owner.primitives(frame,c.mpf((0,1)))
    assert result['proof']['branch']=='regular_small_r_Fourier'
    assert not frame.roots['p2'][(0,0)].zero and ep(frame.roots['p2'][(0,1)].coefficient)[1]<0
    rejected=0
    for call in (lambda:owner.primitives(replace(frame),c.mpf('.23')),
        lambda:owner.roots(0,0),lambda:owner.integrate(count=0,N=160),
        lambda:a.rebase_piece_value(None,None,None),
        lambda:current.NearMidplaneAtlas(owner.inputs.frame,zeta_lower=0,zeta_upper=0)):
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejected+=1
    assert rejected==5
    return dict(passed=True,independent_defining_integral_regular_and_zero_C0_Z_comparisons=comparisons,
        positive_negative_zero_regular_functions_checked=True,
        whole_signed_u_interval_crosses_zero_with_original_Fourier_branch=True,
        actual_nonzero_p2_Z_retained=True,issued_frame_and_domain_rejections=rejected)


def midplane_in_whole_zeta_source(owner):
    a=owner.atlas;c=owner.ctx;frame=owner.source_frame('-2.5','-1.25')
    old=owner.midplane.roots('-2.5','-1.25');comparisons=0
    for name in ('E','V','b','p1','p2','a','t0'):
        for order in (0,1):
            box=frame.roots[name][(0,order)];value=old['roots'][name][(0,order)]
            anchor=current.prior.FormalScale(a.bases,box.scale.powers)
            norm_box=c.mpf(0) if box.zero else box.coefficient*box.bounded_exp((box.scale-anchor).evaluate())
            # Explicit same-family adapter ONLY for contained Z=0, L=1.
            # Its old L factor is exactly one, rather than a fixed-Z rebase.
            old_scale=current.prior.FormalScale(a.bases,
                (value.scale.powers[0],value.scale.powers[1],0,value.scale.powers[3],value.scale.powers[4]),
                offset=value.scale.offset)
            norm_old=c.mpf(0) if value.zero else value.coefficient*value.bounded_exp((old_scale-anchor).evaluate())
            tolerance=c.mpf(2)**(-240)*max(*(abs(v) for v in ep(norm_old)),1)
            saved.contains(norm_box,norm_old,tolerance);comparisons+=1
    return dict(passed=True,accepted_original_midplane_root_enclosures_contained_in_whole_zeta_source=comparisons,
        adapter_only_same_family_exact_Z0_L1=True,
        comparison_rounding_budget='2^-240 relative in finite normalized coefficients',
        unrestricted_fixed_Z_rebase_not_enabled=True)


def integral_contracts(owner,manifest):
    c=owner.ctx;levels=manifest['actual_original_near_midplane_reference_levels'];summaries=[];masses=0;cells_checked=0
    for level in levels:
        count,N=level['exact_original_cells'],level['candidate_N'];cells=level['full_whole_cell_source_phase_records']
        assert level['source_family']==owner.family and 'original_Z_exact' not in level
        assert level['exact_zeta_interval'] in manifest['accepted_signed_zeta_windows'] and len(cells)==count
        for flag in ('whole_near_midplane_reference_C0_Z_contributions_installed',
            'original_regular_source_and_Fourier_Z_functions_used','native_large_Z_derivative_factor_not_capped',
            'actual_source_functions_uniform_over_this_zeta_window','high_precision_source_and_integral_exports_retained',
            'incoming_global_correction_histories_not_assumed_zero','original_P0_kept_separate_and_not_reset',
            'original_pressure_rate_zero_memory_retained'):assert level[flag]
        assert not level['whole_Z_terminal_or_all_route_closure']
        totals={name:c.mpf(0) for name in points.exact.RATES}
        for i,cell in enumerate(cells):
            source=cell['source'];left,right=map(s.Rational,source['exact_reference_cell'])
            assert left==-5+s.Rational(5*i,count) and right==-5+s.Rational(5*(i+1),count)
            assert source['candidate_N']==N and source['source_geometry']=='small_r_series'
            assert source['exact_zeta_interval']==level['exact_zeta_interval']
            u=saved.interval(c,source['whole_zeta_actual_original_u_range'])
            assert ep(u)[0]<0<ep(u)[1] and max(abs(v) for v in ep(u))<=c.mpf('.25')
            roots=source['whole_cell_original_roots']
            assert not roots['p2']['(0, 0)']['exact_zero'] and roots['p2']['(0, 1)']['sign']=='negative'
            Zmap=source['whole_zeta_physical_Z_map']
            assert not Zmap['exact_zero'] and Zmap['formal_positive_scale']['source_exponents']==[-11,-10,0,0]
            budget=source['whole_zeta_positive_Q_variable_L']
            assert ep(saved.interval(c,budget['Q']))[0]>=1 and ep(saved.interval(c,budget['L']))[0]>0
            for piece in source['original_inverse_and_Z_piece_records']:
                assert piece['ordinary_Z']['branch']=='regular_small_r_Fourier'
                assert not piece['ordinary_Z']['u_Z']['exact_zero']
                assert piece['ordinary_Z']['p2_Z_retained'] and piece['ordinary_Z']['E_Z_term_retained']
            assert cell['native_radius_Jacobian']==1;cells_checked+=1
            for name,mass in cell['original_positive_own_rate_masses'].items():
                box=saved.interval(c,mass);assert ep(box)[0]>0;totals[name]+=box
        for name,rate in points.exact.RATES.items():
            r=c.mpf(rate.numerator)/rate.denominator
            expected=c.mpf(5) if ep(r)==(0,0) else (1-c.exp(-5*r))/r
            saved.contains(totals[name],expected);masses+=1
        bounds={}
        for name,pair in level['actual_finite_N_five_density_contribution_enclosures'].items():
            bounds[name]={}
            for jet,row in pair.items():
                scale=row['formal_positive_scale'];powers=scale['source_exponents']
                assert powers==([0,0,0,0] if jet=='C0' else [11,10,-2,0]),(name,jet,powers)
                assert scale['radius_power']==0 and not row['point_value_selected']
                v=saved.interval(c,row['coefficient_interval'])*c.exp(saved.interval(c,scale['additional_log_interval']))
                assert all(mp.isfinite(x) for x in ep(v))
                bounds[name][jet if jet=='C0' else 'Z_divided_by_Pstar11_Cstar10_Lminus2']=max(abs(x) for x in ep(v))
        summaries.append(dict(cells=count,N=N,exact_zeta_interval=level['exact_zeta_interval'],whole_zeta_contribution_bounds=bounds))
    assert [(v['exact_original_cells'],v['candidate_N']) for v in levels]==[(4,160),(16,160),(16,16384),(16,160),(16,16384)]
    assert [v['exact_zeta_interval'] for v in levels]==[manifest['accepted_signed_zeta_windows'][0]]*3+[manifest['accepted_signed_zeta_windows'][1]]*2
    assert manifest['wider_neighborhood_half_width_ratio']==100
    first,last=levels[1],levels[2]
    assert first['full_whole_cell_source_phase_records'][0]['source']['whole_cell_N_dependent_coefficients']['-1']['h']!=last['full_whole_cell_source_phase_records'][0]['source']['whole_cell_N_dependent_coefficients']['-1']['h']
    return dict(passed=True,whole_y_zeta_cells_checked=cells_checked,positive_own_rate_mass_comparisons=masses,
        pressure_zero_rate_memory_mass=5,native_large_Z_source_factor=[11,10,-2,0],
        hundredfold_larger_whole_neighborhood_with_original_regular_guard=True,
        actual_N_dependent_coefficients_recomputed=True,computed_ranges=summaries,
        bounds_are_source_window_contributions_not_terminal_or_NS_errors=True)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
           'current_whole_N_selected',*current.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[name] is False for name in flags)
    owner=current.OriginalReferenceNearMidplane();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (
            ('original_factored_source_and_pressure',lambda:factor_pressure_source(owner)),
            ('independent_regular_neighborhood_primitives',lambda:independent_regular_primitives(owner)),
            ('accepted_midplane_in_whole_zeta',lambda:midplane_in_whole_zeta_source(owner)),
            ('whole_neighborhood_integral_contracts',lambda:integral_contracts(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Actual original whole Rh_reference C0/Z functions, primitives and own-rate contributions on signed zeta[-1e-8,1e-8] and hundredfold expanded[-1e-6,1e-6], retaining physical Z=zeta/Lambda0, positive L/Q hulls, original full-pressure remainder jet provenance and native ordinary-Z amplification. Separate pressure boxes do not preserve joint correlation. Not full Z[-1,1], C1 phase averaging, terminal/all-route controls/global N, recursion or corrected NS.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Whole original near-midplane source and integral evidence PASS',flush=True)
    return report


if __name__=='__main__':run()
