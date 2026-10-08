"""Original regular mixed source, defining integrals and C1 contribution evidence."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_near_midplane_mixed_C1_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;C0,Y,Z,YZ=current.mixed.ORDERS


def original_mixed_source(owner):
    original=owner.owner.templates;inputs=original['inputs'];y=s.Symbol('y',real=True);f=s.exp(y/10)
    sub=dict(zip(inputs[1:5],(f,s.Rational(5,8)*f,s.Rational(5,12)*f*f,s.Rational(5,2)*f*f),strict=True))
    identities=0
    for key,rows in owner.templates.items():
        for (powers,expression),(dpowers,derivative) in zip(original['rows'][key],rows,strict=True):
            assert powers==dpowers
            expected=s.exp(powers[0]*y)*derivative.subs(sub,simultaneous=True)
            function=s.exp(powers[0]*y)*expression.subs(sub,simultaneous=True)
            assert s.simplify(s.diff(function,y)-expected)==0
            for P in original['pressure_symbols']:
                for Q in original['pressure_symbols']:assert s.diff(derivative,P,Q)==0
            identities+=1
    p=mp.mp.clone();p.dps=130;Pstar=p.exp(3);delta=p.exp(-42);Lambda=p.exp(53)
    alpha=p.mpf('2.9')-p.exp(p.mpf('.6'))*p.mpf('.3')/Pstar
    funcs={key:tuple(s.lambdify(inputs,expr,modules=[{'mpf':p.mpf},'mpmath'])
        for powers,expr in rows) for key,rows in original['rows'].items()}
    dfns={key:tuple(s.lambdify(inputs,expr,modules=[{'mpf':p.mpf},'mpmath'])
        for powers,expr in rows) for key,rows in owner.templates.items()}
    def values(y,z):
        f=p.exp(y/10);q=1+z*z
        return (z,f,p.mpf(5)/8*f,p.mpf(5)/12*f*f,p.mpf(5)/2*f*f,
            -alpha/q**2,4*alpha*z/q**3,alpha*(4-20*z*z)/q**4)
    def root(name,y,z,order,radial=False):
        R=110*p.exp(50+y);L=1-delta*z*z
        rows=(owner.templates if radial else original['rows'])[(name,order)]
        fns=(dfns if radial else funcs)[(name,order)]
        return sum(R**r*Pstar**ps*delta**d*L**ell*fn(*values(y,z))
            for ((r,ps,d,ell),expr),fn in zip(rows,fns,strict=True))
    comparisons=0
    for yy in ('-2.337','-.513'):
        yy=p.mpf(yy)
        for zz in ('-.7','0','.7'):
            z=p.mpf(zz)/Lambda
            for name in ('E','V','b','p1','p2'):
                expectedY=p.diff(lambda y:root(name,y,z,0),yy)
                expectedYZ=p.diff(lambda y,z:root(name,y,z,0),(yy,z),(1,1))
                for order,expected in ((0,expectedY),(1,expectedYZ)):
                    got=root(name,yy,z,order,True)
                    assert abs(got-expected)<=p.mpf('1e-105')*max(abs(got),abs(expected),1),(name,yy,zz,order)
                    comparisons+=1
    z,P0Z=inputs[0],original['pressure_symbols'][1];zero={z:0,P0Z:0}
    assert all(s.simplify(expr.subs(zero,simultaneous=True))==0 for powers,expr in owner.templates[('p2',0)])
    assert any(s.simplify(expr.subs(zero,simultaneous=True))!=0 for powers,expr in owner.templates[('p2',1)])
    frame=owner.source_frame('-2.5','-1.25')
    assert not frame.roots['p2'][Y].zero and ep(frame.roots['p2'][YZ].coefficient)[1]<0
    assert all(frame.roots['a'][key].zero and frame.roots['t0'][key].zero for key in (Y,Z,YZ))
    rejected=0
    for call in (lambda:owner.primitive(replace(frame),owner.ctx.mpf('.3')),
        lambda:owner.source_frame(0,0),lambda:owner.integrate(count=0,N=160)):
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejected+=1
    assert rejected==3
    return dict(passed=True,exact_original_C0_Z_row_radial_derivative_identities=identities,
        independent_original_C0_y_and_yZ_derivative_comparisons=comparisons,
        consistent_nonzero_even_pressure_remainder_in_diagnostic=True,
        exact_midplane_u_y_zero_and_u_yZ_nonzero_source_identity=True,
        whole_zeta_nonzero_p2_y_and_negative_p2_yZ_retained=True,
        issued_frame_and_domain_rejections=rejected,finite_units_not_installed_as_original_parameters=True)


def symbolic_regular_and_tail_proof(owner):
    y,z,r,psi=s.symbols('y z r psi',real=True)
    u=s.Function('u')(y,z);h=(1+u*u)**s.Rational(-1,2)
    uy,uz,uyz=s.diff(u,y),s.diff(u,z),s.diff(u,y,z)
    identities=(s.diff(u*h,y,z)-h**3*uyz+3*u*h**5*uy*uz,
        s.diff(h,y,z)+u*h**3*uyz-(3*u*u*h**5-h**3)*uy*uz)
    assert all(s.simplify(v)==0 for v in identities)
    k=s.Symbol('k',integer=True,positive=True)
    # Factor out r^(k-2) for integer k>=2; the remaining polynomial
    # identity extends through r=0 without a negative power evaluation.
    assert s.expand(4*r*r+2*(k-1)*(1-r*r)-
        (2*(k-1)+(6-2*k)*r*r))==0
    t,ty,tz,tp,T2y,T2z,T2yz,T1yz=s.symbols('t ty tz tp T2y T2z T2yz T1yz')
    D=1+t*t;py=-T2y/D;pz=-T2z/D
    pyz=-(T2yz+2*t*ty*pz+2*t*tz*py+2*t*tp*py*pz)/D
    combined=-T2yz/D+2*t*(ty*T2z+tz*T2y)/D**2-2*t*tp*T2y*T2z/D**3
    assert s.cancel(pyz-combined)==0
    total=T1yz+ty*pz+tz*py+tp*py*pz+t*pyz
    expected=T1yz-t*T2yz/D+(t*t-1)*(ty*T2z+tz*T2y)/D**2+(1-t*t)*tp*T2y*T2z/D**3
    assert s.cancel(total-expected)==0
    # Closed geometric tail sums, proved by differentiation.
    R=s.Symbol('R',positive=True);n=s.Symbol('n',integer=True,nonnegative=True)
    G=R**n/(1-R)
    assert s.simplify(R*s.diff(G,R)+G-
        R**n*((n+1)/(1-R)+R/(1-R)**2))==0
    quadratic=R**2*s.diff(G,R,2)+4*R*s.diff(G,R)+2*G
    assert s.simplify(quadratic-R**n*((n+1)*(n+2)/(1-R)+
        (2*n+3)*R/(1-R)**2+R*(1+R)/(1-R)**3))==0
    c=owner.ctx;bounds=current.regular_fourier_tail_bounds(c,c.mpf('.25'))
    assert all(ep(v)[0]>=0 for rows in bounds.values() for v in rows)
    zero=current.regular_fourier_tail_bounds(c,0)
    assert all(ep(v)==(0,0) for rows in zero.values() for v in rows)
    # Bind the constant-nu implicit identity to the ACTUAL reference q
    # definition, not just to an untyped constant-q assumption. Delta=a-2
    # is negative, so the original cutoff is exactly one. The selected
    # eta/log-a datum values are original fixed family parameters.
    frame=owner.source_frame('-2.5','-1.25');aa=frame.roots['a'][C0]
    anchor_checks=0
    for key in (Y,Z,YZ):
        kernel=frame.query['kernel']
        raw=(frame.roots['p2'][key]*kernel.q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate())
        anchored=current.finite_offset_anchor(owner.atlas,raw)
        assert anchored.scale.powers==raw.scale.powers and ep(anchored.scale.offset)==(0,0)
        saved.contains(anchored.coefficient,raw.coefficient*raw.bounded_exp(raw.scale.offset))
        anchor_checks+=1
    delta=aa-2
    assert ep(delta.coefficient)[1]<0
    qa=current.base.current.q_enclosure(aa,delta,
        owner.atlas.copy_interval(owner.owner.scales.logs['eta']),
        owner.atlas.copy_interval(owner.owner.scales.logs['a_min']))
    assert qa['branch']=='active' and ep(qa['sigma_interval'])==(1,1)
    qdiff=owner.atlas.add(qa['q'],-frame.query['kernel'].q)
    assert qdiff.zero or ep(qdiff.coefficient)[0]<=0<=ep(qdiff.coefficient)[1]
    eta=s.Symbol('eta',positive=True)
    qdef=s.sqrt((2*eta-(s.Rational(4,5)-2))/(2*s.Rational(4,5)))
    assert all(s.diff(qdef,y,j,z,k)==0 for j,k in (Y,Z,YZ))
    assert all(s.diff(1+2*qdef*qdef,y,j,z,k)==0 for j,k in (Y,Z,YZ))
    # Tail majorants dominate the absolute derivative coefficient sums for
    # every omitted k>M. The signed polynomial coefficient/k magnitudes
    # are <=2; derivatives use the displayed geometric S1/S2 identities.
    kk=49+n
    gaps=(1-(kk-1)/kk,(kk-2)-(kk-1)*(kk-2)/kk,
        2-2*(kk-1)/kk,2-(2*kk-6)/kk,
        2*(kk-2)-2*(kk-1)*(kk-2)/kk,
        2*(kk-2)*(kk-3)-2*(kk-1)*(kk-2)*(kk-3)/kk,
        2*kk-(2*kk-6),2*kk*(kk-1)-(2*kk-6)*(kk-1))
    assert all(s.factor(gap).is_positive is True for gap in gaps)
    return dict(passed=True,smooth_u_mixed_chain_identities=2,original_T2_polynomial_identity=True,
        true_phase_implicit_inverse_and_total_T1_mixed_identities=2,
        exact_first_second_geometric_tail_sum_identities=2,
        tail_majorants='For all integer k>M>=4: (k-1)/k<=1; (k-1)(k-2)/k<=k-2; 2(k-1)/k<=2; (2k-6)/k<=2. Differentiate each original monomial before the bound and sum S1/S2. M48 only nonnegative powers.',
        actual_analytic_tail_partials_not_derivatives_of_radius_caps=True,
        tails_exact_zero_at_r_zero=True,original_reference_q_and_nu_constant_scope=True,
        finite_offset_anchor_source_enclosure_checks=anchor_checks,
        actual_q_definition_negative_Delta_exact_sigma_one_and_nu_slow_jets_bound=True)


def independent_regular_mixed_integrals(owner):
    a=owner.atlas;c=owner.ctx;p=mp.mp.clone();p.dps=110
    theta,uu,qq=s.symbols('theta u q',real=True)
    h=s.sqrt(1+uu*uu);r=uu/h
    direction=2*qq/h*(s.cos(theta)-r)/(1-2*r*s.cos(theta)+r*r)
    integrands={power:tuple(s.lambdify((theta,uu,qq),s.diff(direction**power,uu,order),
        modules=[{'mpf':p.mpf},'mpmath']) for order in (0,1,2)) for power in (1,2)}
    comparisons=0;traces=0
    for utext in ('-.14','0','.14'):
        u0=p.mpf(utext);q=p.mpf('.9');uy=p.mpf('.3');uz=p.mpf('-.4');uyz=p.mpf('.11')
        rows={name:current.MixedJet(a,{key:a.scalar(c.mpf(value)) for key,value in zip(current.mixed.ORDERS,vals,strict=True)})
            for name,vals in dict(a=('.8',0,0,0),b=(0,0,0,0),t0=(0,0,0,0),
                E=('1.3','-.2','.17','.08'),p2=(u0/q,uy/q,uz/q,uyz/q)).items()}
        directional={name:{C0:row[C0],Z:row[Z]} for name,row in rows.items()}
        kernel=current.whole.conditioned.PositiveLogQPhase(
            dict(q=a.scalar(c.mpf(q)),roots=directional,original_u_source=rows['p2'][C0]*a.scalar(c.mpf(q))),c.mpf(0))
        assert kernel.geometry=='small_r_series'
        for psi_text in ('.47','3.9'):
            ps=p.mpf(psi_text);segments=[0]+([p.pi] if ps>p.pi else [])+[ps]
            got=current.regular_fixed_and_implicit_mixed(a,kernel,rows,c.mpf(ps/(2*p.pi)))
            partials={}
            for power,name in ((1,'T1'),(2,'T2')):
                f0,f1,f2=integrands[power]
                partials[name]={C0:p.quad(lambda th:f0(th,u0,q),segments),
                    Y:p.quad(lambda th:f1(th,u0,q)*uy,segments),
                    Z:p.quad(lambda th:f1(th,u0,q)*uz,segments),
                    YZ:p.quad(lambda th:f2(th,u0,q)*uy*uz+f1(th,u0,q)*uyz,segments)}
                record=got['record']['original_fixed_angle_'+name]
                for key,value in partials[name].items():
                    saved.contains(saved.saved_value(c,record[str(key)]),value,p.mpf('1e-88'));comparisons+=1
            tt=integrands[1][0](ps,u0,q);du=integrands[1][1](ps,u0,q)
            ty,tz=du*uy,du*uz
            tp=p.diff(lambda th:integrands[1][0](th,u0,q),ps);D=1+tt*tt
            py=-partials['T2'][Y]/D;pz=-partials['T2'][Z]/D
            pyz=-(partials['T2'][YZ]+2*tt*ty*pz+2*tt*tz*py+2*tt*tp*py*pz)/D
            totalY=partials['T1'][Y]+tt*py;totalZ=partials['T1'][Z]+tt*pz
            totalYZ=partials['T1'][YZ]+ty*pz+tz*py+tp*py*pz+tt*pyz
            phi=(ps+partials['T2'][C0])/(2*p.pi*(1+2*q*q));aa=p.mpf('.8')
            expectedA={C0:aa/2*(phi-ps/(2*p.pi)),Y:-aa*py/(4*p.pi),
                Z:-aa*pz/(4*p.pi),YZ:-aa*pyz/(4*p.pi)}
            expectedB={C0:-aa*p.mpf('1.3')*partials['T1'][C0]/(4*p.pi),
                Y:-aa/(4*p.pi)*(p.mpf('-.2')*partials['T1'][C0]+p.mpf('1.3')*totalY),
                Z:-aa/(4*p.pi)*(p.mpf('.17')*partials['T1'][C0]+p.mpf('1.3')*totalZ),
                YZ:-aa/(4*p.pi)*(p.mpf('.08')*partials['T1'][C0]+p.mpf('-.2')*totalZ+
                    p.mpf('.17')*totalY+p.mpf('1.3')*totalYZ)}
            for name,expected in (('A',expectedA),('B',expectedB)):
                for key,value in expected.items():
                    saved.contains(got[name][key].finite_interval(),value,p.mpf('1e-88'));comparisons+=1
        for x in (0,'.5',1):
            got=current.regular_fixed_and_implicit_mixed(a,kernel,rows,c.mpf(x))
            assert all(got[name][key].zero for name in ('A','B') for key in current.mixed.ORDERS);traces+=1
    frame=owner.source_frame('-2.5','-1.25');whole=owner.primitive(frame,c.mpf((0,1)))
    assert whole['record']['no_division_by_r_or_p2']
    assert not frame.query['kernel'].u.zero
    nonconstant=dict(frame.roots);v=nonconstant['a']
    nonconstant['a']=current.MixedJet(a,{C0:v[C0],Y:a.scalar(1),Z:a.scalar(0),YZ:a.scalar(0)})
    try:current.regular_fixed_and_implicit_mixed(a,frame.query['kernel'],nonconstant,c.mpf('.3'))
    except ValueError:pass
    else:raise AssertionError('O2 nonconstant parameter guard missing')
    return dict(passed=True,independent_original_defining_integral_and_primitive_C0_y_Z_yZ_comparisons=comparisons,
        positive_negative_zero_u_checked=True,exact_all_jet_period_halfperiod_traces=traces,
        whole_zeta_crossing_not_declared_exact_midplane=True,
        complete_E_yZ_and_both_first_cross_terms_tested=True,
        O2_nonconstant_reference_parameter_rejected=True)


def normalized_bound(owner,row,jet):
    c=owner.ctx;scale=row['formal_positive_scale'];powers=scale['source_exponents']
    assert powers[:2]==([0,0] if jet==C0 else [11,10]),(jet,powers)
    assert powers[3]==0 and scale['radius_power']==0
    # Only cancel Lambda0; all original variable-L powers remain in a
    # directed bounded normalization, never astronomical native exp.
    log=saved.interval(c,scale['additional_log_interval'])+owner.atlas.bases[2]*powers[2]
    value=saved.interval(c,row['coefficient_interval'])*c.exp(log)
    return max(abs(v) for v in ep(value))


def integral_contracts(owner,manifest):
    c=owner.ctx;levels=manifest['actual_original_near_midplane_mixed_C1_levels'];summaries=[];cells_checked=0;masses=0
    for level in levels:
        count,N=level['exact_radial_cells'],level['candidate_N'];cells=level['whole_cell_mixed_source_and_IBP_records']
        assert len(cells)==count and 'original_Z_exact' not in level
        assert level['source_family']==owner.family and level['source_window']==['-5','0']
        assert level['exact_zeta_interval']==['-1/1000000','1/1000000']
        assert level['actual_mixed_source_and_C1_averaging_installed_on_whole_near_midplane_window']
        assert not level['O2_nonconstant_parameter_extension_installed'] and not level['whole_Z_terminal_or_all_route_closure']
        assert level['actual_global_C0_Z_endpoint_terms_retained'] and level['same_original_source_internal_C0_Z_traces_cancel']
        assert level['no_cross_chart_seam_assumed'] and level['incoming_histories_and_P0_not_reset']
        assert level['pressure_zero_rate_memory_retained'] and level['original_interval_export_work_precision_retained']
        totals={name:c.mpf(0) for name in current.points.exact.RATES}
        for i,cell in enumerate(cells):
            source=cell['source'];frame=source['source'];left,right=map(s.Rational,frame['exact_reference_cell'])
            assert left==-5+s.Rational(5*i,count) and right==-5+s.Rational(5*(i+1),count)
            assert frame['source_family']==owner.family and source['candidate_N']==N
            assert frame['original_L_Z_not_applied_twice'] and frame['original_pressure_datum_y_and_yZ_exact_zero']
            assert source['complete_Z_product_rules_and_exact_N_remainder_retained']
            primitive=source['actual_mixed_primitive']
            for flag in ('no_division_by_r_or_p2','both_fixed_angle_cross_terms_and_implicit_curvature_retained',
                'original_reference_a_b_t0_q_nu_slow_derivatives_exact_zero',
                'finite_radial_source_offset_hulls_enclosed_once_before_Fourier_products',
                'finite_offset_anchor_retains_all_original_parameter_powers',
                'tail_derivatives_are_actual_analytic_tail_partials_not_bound_derivatives'):assert primitive[flag]
            radius=mp.mp.make_mpf(tuple(primitive['whole_regular_r_magnitude_upper']['exact_mpf_tuple']))
            assert primitive['M']==48 and radius<=c.mpf('.25')
            assert primitive['denominator_one_plus_t_squared_lower']==1
            assert not primitive['symmetry_trace_exact_zero']
            assert frame['actual_root_jets']['p2'][str(YZ)]['sign']=='negative'
            assert not frame['actual_root_jets']['p2'][str(Z)]['exact_zero']
            for endpoint in source['actual_original_endpoint_phases'].values():
                assert endpoint['positive_original_origin_offset']['strictly_positive']
            assert source['native_log_radius_Jacobian']==1 and source['phase_Z_exact_zero']
            for name,row in cell['contributions'].items():
                assert row['N_power']==-2
                mass=saved.interval(c,row['positive_own_rate_mass']);assert ep(mass)[0]>0;totals[name]+=mass
                for jet,part in row['C0_Z_rows'].items():
                    if 0<i<count-1:assert part['global_endpoint_term']['exact_zero']
            cells_checked+=1
        for name,rate in current.points.exact.RATES.items():
            rate=c.mpf(rate.numerator)/rate.denominator
            expected=c.mpf(5) if ep(rate)==(0,0) else (1-c.exp(-5*rate))/rate
            saved.contains(totals[name],expected);masses+=1
        bounds={}
        for jet in (C0,Z):
            bounds[str(jet)]={}
            for name in current.points.exact.RATES:
                values={kind:normalized_bound(owner,level['C0_Z_'+kind+'_contribution_enclosures'][str(jet)][name],jet)
                    for kind in ('phase_averaged','direct','effective')}
                assert values['effective']<=min(values['phase_averaged'],values['direct'])
                bounds[str(jet)][name]=values
        summaries.append(dict(cells=count,N=N,C0_bounds_and_Z_bounds_divided_by_Lambda0=bounds,
            selections=level['effective_bound_selection']))
    assert [(v['exact_radial_cells'],v['candidate_N']) for v in levels]==[(4,160),(16,160),(16,16384)]
    first,last=levels[1],levels[-1]
    assert first['whole_cell_mixed_source_and_IBP_records'][0]['source']['actual_N_nonlinear_C0_Z_remainder']!=last['whole_cell_mixed_source_and_IBP_records'][0]['source']['actual_N_nonlinear_C0_Z_remainder']
    return dict(passed=True,actual_whole_regular_mixed_cells_checked=cells_checked,
        positive_original_own_rate_mass_comparisons=masses,pressure_zero_rate_memory_mass=5,
        actual_N_dependent_nonlinear_coefficients_recomputed=True,complete_C0_Z_IBP_and_direct_selection=True,
        normalized_Z_divisor='Lambda0=Pstar^11*Cstar^10; actual L powers retained',
        computed_ranges=summaries,terminal_or_global_N_not_claimed=True)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*current.near.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalReferenceNearMidplaneMixed();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (
            ('original_radial_and_mixed_source',lambda:original_mixed_source(owner)),
            ('regular_polynomial_chain_and_tail_proof',lambda:symbolic_regular_and_tail_proof(owner)),
            ('independent_regular_mixed_defining_integrals',lambda:independent_regular_mixed_integrals(owner)),
            ('actual_whole_neighborhood_C1_contributions',lambda:integral_contracts(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Original reference whole signed zeta[-1e-6,1e-6] actual C0/y/Z/yZ source and regular Fourier/implicit mixed primitive functions, complete E/V products and actual-N C0/Z own-rate IBP/direct bounds. Constant reference q/nu only. Native Z factor remains, no full-Z/O2/all-route/control/global frequency/recursion/corrected-NS closure.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original whole near-midplane mixed source and C1 contribution evidence PASS',flush=True)
    return report


if __name__=='__main__':run()
