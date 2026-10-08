"""Independent mixed source, implicit-chain and weighted C1 integration evidence."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_mixed_C1_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep
C0,Y,Z,YZ=current.ORDERS


def original_mixed_source_references(owner):
    inputs=owner.parent.owner.templates['inputs'];y=s.Symbol('y',real=True);f=s.exp(y/10)
    substitutions=dict(zip(inputs[1:5],(f,s.Rational(5,8)*f,s.Rational(5,12)*f*f,s.Rational(5,2)*f*f),strict=True))
    identities=0
    for name,rows in owner.templates.items():
        ordinary=owner.parent.owner.templates['rows'][(name,1)]
        for (powers,expr),(mixed_powers,mixed) in zip(ordinary,rows,strict=True):
            assert powers==mixed_powers
            original=s.exp(powers[0]*y)*expr.subs(substitutions,simultaneous=True)
            expected=s.exp(powers[0]*y)*mixed.subs(substitutions,simultaneous=True)
            assert s.simplify(s.diff(original,y)-expected)==0,(name,powers)
            for P in owner.parent.owner.templates['pressure_symbols']:
                for Q in owner.parent.owner.templates['pressure_symbols']:
                    assert s.diff(mixed,P,Q)==0
            identities+=1
    # Differentiate the original C0 function in BOTH variables, including
    # the original L(Z) and shared analytic P0(Z). No Z-row reuse here.
    p=mp.mp.clone();p.dps=140
    functions={}
    for (name,order),rows in owner.parent.owner.templates['rows'].items():
        if order==0:
            functions[name]=tuple(s.lambdify(inputs,expr,modules=[{'mpf':p.mpf},'mpmath']) for powers,expr in rows)
    mixed_functions={name:tuple(s.lambdify(inputs,expr,modules=[{'mpf':p.mpf},'mpmath']) for powers,expr in rows)
                     for name,rows in owner.templates.items()}
    alpha=p.mpf('2.9');Pstar=p.exp(3);delta=p.exp(-42)
    def values(yy,zz):
        ff=p.exp(yy/10);q=1+zz*zz
        return (zz,ff,p.mpf(5)/8*ff,p.mpf(5)/12*ff*ff,p.mpf(5)/2*ff*ff,
                -alpha/q**2,4*alpha*zz/q**3,alpha*(4-20*zz*zz)/q**4)
    def root(name,yy,zz,mixed=False):
        R=110*p.exp(50+yy);L=1-delta*zz*zz
        rows=owner.templates[name] if mixed else owner.parent.owner.templates['rows'][(name,0)]
        fns=mixed_functions[name] if mixed else functions[name]
        return sum(R**r*Pstar**a*delta**d*L**ell*fn(*values(yy,zz))
                   for ((r,a,d,ell),expr),fn in zip(rows,fns,strict=True))
    comparisons=0
    for yy,zz in (('-2.337','.37'),('-.513','-.41')):
        yy,zz=p.mpf(yy),p.mpf(zz)
        for name in ('E','V','b','p1','p2'):
            expected=p.diff(lambda y,z:root(name,y,z),(yy,zz),(1,1))
            got=root(name,yy,zz,True)
            assert abs(expected-got)<=p.mpf('1e-110')*max(abs(expected),abs(got),1),(name,yy,zz)
            comparisons+=1
    frame=owner.source_frame('-2.5','-1.25')
    for name in ('V','b','a','t0'):assert frame.roots[name][YZ].zero
    assert not frame.roots['E'][YZ].zero and not frame.roots['p2'][YZ].zero
    rejected=0
    for call in (lambda:owner.primitive(replace(frame),owner.ctx.mpf('.3')),
                 lambda:owner.source_frame(0,0),
                 lambda:current.OriginalReferenceMixedC1(Z=0),
                 lambda:owner.integrate(count=0,N=160)):
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejected+=1
    assert rejected==4
    return dict(passed=True,exact_original_Z_row_radial_derivative_identities=identities,
        independent_original_C0_yZ_function_derivative_comparisons=comparisons,
        original_L_Z_applied_once_and_shared_pressure_Z_jets_retained=True,
        pressure_affine_late_factor_retained=True,issued_frame_and_domain_rejections=rejected,
        finite_diagnostic_units_not_installed_in_original_frame=True)


def symbolic_chain_and_ranges():
    y,z=s.symbols('y z',real=True)
    u=s.Function('u')(y,z);h=(1+u*u)**s.Rational(-1,2);r=u*h;ss=h*h
    uy,uz,uyz=s.diff(u,y),s.diff(u,z),s.diff(u,y,z)
    KY,KZ,KYZ=uy*h,uz*h,uyz*h
    identities=(
        s.diff(r,y,z)-ss*(KYZ-3*r*KY*KZ),
        s.diff(ss,y,z)+2*ss*(r*KYZ+(ss-3*r*r)*KY*KZ),
        s.diff(h,y,z)-h*((3*r*r-1)*KY*KZ-r*KYZ))
    assert all(s.simplify(expr)==0 for expr in identities)
    # Independent substitution in the defining implicit mixed equations.
    t,ty,tz,tp,T2y,T2z,T2yz,T1yz=s.symbols('t ty tz tp T2y T2z T2yz T1yz')
    D=1+t*t;py=-T2y/D;pz=-T2z/D
    pyz=-(T2yz+2*t*ty*pz+2*t*tz*py+2*t*tp*py*pz)/D
    expected=-T2yz/D+2*t*(ty*T2z+tz*T2y)/D**2-2*t*tp*T2y*T2z/D**3
    assert s.cancel(pyz-expected)==0
    total=T1yz+ty*pz+tz*py+tp*py*pz+t*pyz
    combined=T1yz-t*T2yz/D+(t*t-1)*(ty*T2z+tz*T2y)/D**2+(1-t*t)*tp*T2y*T2z/D**3
    assert s.cancel(total-combined)==0
    by,bz,gy,gz=s.symbols('by bz gy gz')
    assert s.expand((ty*T2z+tz*T2y).subs({ty:by*t+gy,tz:bz*t+gz})-
                    (t*(by*T2z+bz*T2y)+gy*T2z+gz*T2y))==0
    # Exact nonnegative square certificates for every nontrivial rational range.
    t=s.Symbol('t',real=True);D=1+t*t
    certificates=(
        s.Rational(1,2)-t/D-(t-1)**2/(2*D),
        s.Rational(1,2)+t/D-(t+1)**2/(2*D),
        s.Rational(1,2)-2*t*t/D**2-(t*t-1)**2/(2*D**2),
        (t*t-1)/D**2+1-(t**4+3*t*t)/D**2,
        s.Rational(1,8)-(t*t-1)/D**2-(t*t-3)**2/(8*D**2),
        s.Rational(1,4)-t*(t*t-1)/D**2-(t*t-2*t-1)**2/(4*D**2),
        s.Rational(1,4)+t*(t*t-1)/D**2-(t*t+2*t-1)**2/(4*D**2))
    assert all(s.cancel(expr)==0 for expr in certificates)
    k,n,r,h=s.symbols('k n r h',real=True)
    assert s.expand(5*(1+k*k)-(1+2*k)**2-(k-2)**2)==0
    assert s.expand(20*(1+k*k)-(4+2*k)**2-4*(2*k-1)**2)==0
    assert s.expand(1-(n-r)**2-(1-r*r+2*r*n-n*n))==0
    assert s.simplify((h**s.Rational(3,2)*(1+(n/h)**2)**s.Rational(1,4))**4-
                      h**4*(h*h+n*n))==0
    # Exact N remainder and its ordinary Z derivative.
    x=s.Symbol('x',real=True);R2=(s.exp(x)-1-x)/x**2
    assert s.simplify(2*R2+x*s.diff(R2,x)-(s.exp(x)-1)/x)==0
    return dict(passed=True,exact_r_s_hinv_mixed_chain_identities=3,
        original_implicit_inverse_and_total_T1_mixed_identities=2,
        exact_rational_nonnegative_certificates=len(certificates),
        source_geometry_and_Cauchy_weighted_curvature_certificates=4,
        exact_Q_Z_exprel_identity=True,
        weighted_curvature_proof='With h=hinv, n=r+coschi, k=|n|/h and D=1+t^2>=1+k^2: P=|sinchi|*|Dchi-3*h^2|+h^2*2*pi/|r|. P^2<=2*h^3*(h+2*k)*(4*h+2*k)^2+2*h^4*D0^2 and |sinchi|*Dchi<=h^(3/2)*(h+2*k)^(3/2). Multiply before bounding t_psi*J^2. The displayed square certificates and h^4*(h^2+n^2)<=5 give C25,C2 for denominator powers5/2,2. The rational prefactors then give the two stated source-correlated bounds.')


def independent_mixed_integral_references(owner):
    a=owner.atlas;c=owner.ctx;p=mp.mp.clone();p.dps=120
    th,uu,qq=s.symbols('theta u q',real=True)
    hh=s.sqrt(1+uu*uu);rr=uu/hh
    original=2*qq/hh*(s.cos(th)-rr)/(1-2*rr*s.cos(th)+rr*rr)
    integrands={power:tuple(s.lambdify((th,uu,qq),s.diff(original**power,uu,order),
                 modules=[{'mpf':p.mpf},'mpmath']) for order in (0,1,2)) for power in (1,2)}
    comparisons=0
    for sign in (1,-1):
        q=p.mpf('.9');u0=p.mpf('1.8')*sign;uy=p.mpf('.3');uz=p.mpf('-.4');uyz=p.mpf('.11')
        def direction(theta,y,z):
            uu=u0+uy*y+uz*z+uyz*y*z;hh=p.sqrt(1+uu*uu);rr=uu/hh
            return 2*q/hh*(p.cos(theta)-rr)/(1-2*rr*p.cos(theta)+rr*rr)
        for psi_text in ('.47','3.9'):
            psi=p.mpf(psi_text);rr=u0/p.sqrt(1+u0*u0)
            chi=2*p.atan2((1+rr)*p.sin(psi/2),(1-rr)*p.cos(psi/2))
            if chi<0:chi+=2*p.pi
            coordinate=c.mpf(chi/(2*p.pi))
            rows={}
            for name,vals in dict(a=('.8',0,0,0),b=(0,0,0,0),t0=(0,0,0,0),
                                  E=('1.3','-.2','.17','.08'),p2=(u0/q,uy/q,uz/q,uyz/q)).items():
                rows[name]=current.MixedJet(a,{key:a.scalar(c.mpf(value)) for key,value in zip(current.ORDERS,vals,strict=True)})
            directional={name:{(0,0):row[C0],(0,1):row[Z]} for name,row in rows.items()}
            kernel=current.whole.conditioned.PositiveLogQPhase(
                dict(q=a.scalar(c.mpf(q)),roots=directional,original_u_source=rows['p2'][C0]*a.scalar(c.mpf(q))),c.mpf(0))
            got=current.fixed_and_implicit_mixed(a,kernel,rows,coordinate)
            # Independently integrate the ORIGINAL rational direction and its
            # squared function at fixed psi, then apply the defining chain.
            partials={}
            segments=[0]+([p.pi] if psi>p.pi else [])+[psi]
            for power,name in ((1,'T1'),(2,'T2')):
                f0,f1,f2=integrands[power]
                partials[name]={C0:p.quad(lambda theta:f0(theta,u0,q),segments),
                    Y:p.quad(lambda theta:f1(theta,u0,q)*uy,segments),
                    Z:p.quad(lambda theta:f1(theta,u0,q)*uz,segments),
                    YZ:p.quad(lambda theta:f2(theta,u0,q)*uy*uz+f1(theta,u0,q)*uyz,segments)}
                record=got['record']['original_fixed_angle_'+name]
                for key,value in partials[name].items():
                    try:saved.contains(saved.saved_value(c,record[str(key)]),value,p.mpf('1e-88'))
                    except AssertionError:
                        box=saved.saved_value(c,record[str(key)])
                        raise AssertionError((sign,psi_text,name,key,p.nstr(ep(box)[0]-value,30),p.nstr(ep(box)[1]-value,30)))
                    comparisons+=1
            tt=direction(psi,0,0)
            ty=p.diff(lambda y:direction(psi,y,0),0);tz=p.diff(lambda z:direction(psi,0,z),0)
            tp=p.diff(lambda theta:direction(theta,0,0),psi);D=1+tt*tt
            h=1/p.sqrt(1+u0*u0);ss=h*h;n=rr+p.cos(chi);Dchi=ss+2*rr*n
            J=2*q*q/rr**2*(p.sin(chi)*(Dchi-3*ss)+ss*(chi-psi)/rr)
            assert abs(partials['T2'][Y]-uy*h*J)<=p.mpf('1e-88')
            assert abs(partials['T2'][Z]-uz*h*J)<=p.mpf('1e-88')
            assert abs(ty-(uy*h*(2*p.cos(chi)-rr)*tt-2*q*uy*h*h))<=p.mpf('1e-88')
            assert abs(tp+2*q*p.sin(chi)*Dchi/h**3)<=p.mpf('1e-88')
            comparisons+=4
            py=-partials['T2'][Y]/D;pz=-partials['T2'][Z]/D
            pyz=-(partials['T2'][YZ]+2*tt*ty*pz+2*tt*tz*py+2*tt*tp*py*pz)/D
            totalY=partials['T1'][Y]+tt*py;totalZ=partials['T1'][Z]+tt*pz
            totalYZ=partials['T1'][YZ]+ty*pz+tz*py+tp*py*pz+tt*pyz
            AYZ=-p.mpf('.8')*pyz/(4*p.pi)
            BYZ=-p.mpf('.8')/(4*p.pi)*(p.mpf('.08')*partials['T1'][C0]-
                  p.mpf('.2')*totalZ+p.mpf('.17')*totalY+p.mpf('1.3')*totalYZ)
            saved.contains(got['A'][YZ].finite_interval(),AYZ,p.mpf('1e-88'))
            saved.contains(got['B'][YZ].finite_interval(),BYZ,p.mpf('1e-88'));comparisons+=2
    # Extreme finite diagnostic u tests supplement the algebraic universal
    # proof. They are not the proof and never replace original source values.
    curvature_tests=0
    for uu in ('-.01','.01','-1.8','1.8','-20','20','-1e6','1e6'):
        uu=p.mpf(uu);h=1/p.sqrt(1+uu*uu);r=uu*h
        for q in (p.mpf('.5'),p.mpf('.9'),p.mpf(2)):
            D0=2*p.pi/abs(r)
            C25=2*(20*5**(p.mpf(5)/4)+5**(p.mpf(3)/4)*D0**2)
            C2=2*(20*5**(p.mpf(3)/2)+5**(p.mpf(3)/4)*D0**2)
            for fraction in ('.000001','.07','.23','.49','.5','.500001','.71','.999999'):
                chi=2*p.pi*p.mpf(fraction)
                psi=2*p.atan2((1-r)*p.sin(chi/2),(1+r)*p.cos(chi/2))
                if psi<0:psi+=2*p.pi
                n=r+p.cos(chi);Dchi=h*h+2*r*n;t=2*q*n/h
                tp=-2*q*p.sin(chi)*Dchi/h**3
                J=2*q*q/r**2*(p.sin(chi)*(Dchi-3*h*h)+h*h*(chi-psi)/r)
                D=1+t*t
                assert abs(2*t*tp*J*J/D**3)<=16*q**5/r**4*C25
                assert abs((1-t*t)*tp*J*J/D**3)<=8*q**5/r**4*C2
                curvature_tests+=2
    return dict(passed=True,independent_original_defining_integral_mixed_comparisons=comparisons,
        independent_finite_curvature_diagnostics=curvature_tests,positive_and_negative_r=True,
        near_small_r_and_large_u_diagnostics=True,universal_bound_from_algebra_not_samples=True)


def actual_integral_contracts(owner,manifest):
    c=owner.ctx;levels=manifest['actual_original_mixed_C1_reference_levels'];ranges=[]
    for level in levels:
        count,N=level['exact_radial_cells'],level['candidate_N'];cells=level['whole_cell_mixed_source_and_IBP_records']
        assert len(cells)==count and level['source_family']==owner.family
        assert level['actual_global_C0_Z_endpoint_terms_retained']
        assert level['same_original_source_internal_C0_Z_traces_cancel']
        assert level['no_cross_chart_seam_assumed'] and level['incoming_histories_and_P0_not_reset']
        assert level['pressure_zero_rate_memory_retained'] and not level['whole_Z_terminal_or_all_route_closure']
        assert level['original_interval_export_work_precision_retained']
        totals={name:c.mpf(0) for name in current.points.exact.RATES}
        for i,cell in enumerate(cells):
            source=cell['source'];frame=source['source']
            left,right=map(s.Rational,frame['exact_reference_cell'])
            assert left==-5+s.Rational(5*i,count) and right==-5+s.Rational(5*(i+1),count)
            assert frame['source_family']==owner.family and source['candidate_N']==N
            assert frame['original_L_Z_not_applied_twice'] and frame['shared_P0_y_and_yZ_exact_zero']
            assert frame['source_pressure_late_factors_retained']
            assert source['native_log_radius_Jacobian']==1 and source['phase_Z_exact_zero']
            assert source['complete_Z_product_rules_and_exact_N_remainder_retained']
            primitive=source['actual_mixed_primitive']
            assert primitive['both_fixed_angle_first_cross_terms_and_curvature_retained']
            assert primitive['source_enclosures_not_saved_cap_or_midpoint_values']
            curve=primitive['correlated_weighted_curvature']
            assert curve['microscopic_source_factors_collected_before_bounding']
            assert curve['cap_not_selected_as_field_point_value']
            assert ep(saved.interval(c,curve['original_q_directed_range']))[0]>=c.mpf('.5')
            for endpoint in source['actual_original_endpoint_phases'].values():
                assert endpoint['positive_original_origin_offset']['strictly_positive']
            for name,row in cell['contributions'].items():
                assert row['N_power']==-2
                totals[name]+=saved.interval(c,row['positive_own_rate_mass'])
                for jet,record in row['C0_Z_rows'].items():
                    if 0<i<count-1:assert record['global_endpoint_term']['exact_zero']
        saved.contains(totals['p'],c.mpf(5))
        values={}
        for jet in (C0,Z):
            values[str(jet)]={}
            for name in ('h','e','p'):
                pair={typ:max(abs(v) for v in ep(saved.saved_value(c,
                      level['C0_Z_'+typ+'_contribution_enclosures'][str(jet)][name],owner.atlas.bases)))
                      for typ in ('phase_averaged','direct','effective')}
                assert pair['effective']<=min(pair['phase_averaged'],pair['direct'])
                values[str(jet)][name]=pair
        ranges.append(dict(cells=count,N=N,C0_Z_normalized_absolute_bounds=values))
    assert [(v['exact_radial_cells'],v['candidate_N']) for v in levels]==[(4,160),(16,160),(16,320),(16,16384)]
    for row in ranges[-1]['C0_Z_normalized_absolute_bounds'][str(Z)].values():
        assert row['phase_averaged']<row['direct']/20
    assert all(level['effective_bound_selection'][str(Z)]['p']=='direct_smaller_by_collected_source_factor_ratio'
               for level in levels[0:3])
    # Actual N-specific nonlinear source records must change, not be relabelled.
    first,last=levels[1],levels[-1]
    assert first['whole_cell_mixed_source_and_IBP_records'][0]['source']['actual_N_nonlinear_C0_Z_remainder']!=last['whole_cell_mixed_source_and_IBP_records'][0]['source']['actual_N_nonlinear_C0_Z_remainder']
    # The same accepted C0 method is recomputed. Its collected formal anchor
    # can differ, and the historical export used a lower real work precision.
    # Compare physical bounds with a separately declared diagnostic budget;
    # this budget never enters an installed field or integral enclosure.
    old=json.loads(gzip.decompress((current.HERE/current.phase.NAME).read_bytes()))
    budget=c.mpf(2)**-48
    for level,prior in zip(levels[:3],old['actual_original_reference_phase_averaged_levels'],strict=True):
        for name in current.points.exact.RATES:
            new=level['C0_Z_effective_contribution_enclosures'][str(C0)][name]
            previous=prior['C0_effective_contribution_enclosures'][name]
            ns,ps=new['formal_positive_scale'],previous['formal_positive_scale']
            assert ns['source_exponents']==ps['source_exponents'] and ns['radius_power']==ps['radius_power']
            # Cancel the identical original factor before normalization;
            # microscopic m/k bounds must never materialize that exponential.
            new_bound=ep(saved.interval(c,new['coefficient_interval'])*c.exp(
                saved.interval(c,ns['additional_log_interval'])))[1]
            old_bound=ep(saved.interval(c,previous['coefficient_interval'])*c.exp(
                saved.interval(c,ps['additional_log_interval'])))[1]
            assert abs(new_bound-old_bound)<=ep(budget*c.mpf(old_bound))[1]
            # Nonzero exported coefficients must retain original interval
            # precision, not a rounded double-precision endpoint.
            if name in ('h','e','p') and not new['exact_zero']:
                bits=new['coefficient_interval']['upper_exact_mpf_tuple'][3]
                assert bits>200
    return dict(passed=True,actual_levels=[(v['exact_radial_cells'],v['candidate_N']) for v in levels],
        true_mixed_source_and_fresh_N_dependent_Z_remainders=True,
        N16384_Z_averaging_bounds_over20_times_tighter_than_direct=True,
        coarse_N_direct_Z_bounds_retained_when_tighter=True,
        predecessor_C0_physical_bounds_agree_with_separate_2minus48_relative_comparison_budget=True,
        comparison_budget_not_used_as_source_or_integral_value=True,
        high_precision_directed_interval_exports_retained=True,computed_ranges=ranges,
        pressure_zero_rate_memory_mass=5,whole_Z_or_global_N_not_claimed=True)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
           'current_whole_N_selected',*current.whole.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[name] is False for name in flags)
    owner=current.OriginalReferenceMixedC1(Z='.37')
    with mp.workdps(owner.ctx.dps+40):
        checks={}
        for name,call in (
            ('actual_original_mixed_roots',lambda:original_mixed_source_references(owner)),
            ('exact_implicit_and_weighted_curvature_identities',symbolic_chain_and_ranges),
            ('independent_defining_mixed_source_evidence',lambda:independent_mixed_integral_references(owner)),
            ('actual_C1_source_integral_evidence',lambda:actual_integral_contracts(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
        Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Actual original reference C0,y,Z,yZ roots, correlated implicit mixed primitives, and fixed nonzero-Z C0/Z own-rate phase averaging. Algebraic weighted curvature proof plus independent defining integral diagnostics. Not whole-Z/midplane/all-route terminal closure, selected global N, actual recursion or corrected Cartesian NS.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original mixed source, correlated phase C1 integration and independent references PASS',flush=True)
    return report


if __name__=='__main__':run()
