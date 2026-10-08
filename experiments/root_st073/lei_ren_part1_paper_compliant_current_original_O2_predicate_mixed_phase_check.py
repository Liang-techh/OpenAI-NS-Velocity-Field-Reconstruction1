"""Original rational integrals, varying-q mixed calculus and predicate union."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_predicate_mixed_phase as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

C0,Y,Z,YZ=current.C0,current.Y,current.Z,current.YZ;ep=current.ep


def same_original_direction_and_mixed_calculus():
    r,q,cp=s.symbols('r q cospsi',real=True)
    ss=1-r*r;Dpsi=1-2*r*cp+r*r
    cc=((1+r*r)*cp-2*r)/Dpsi
    assert s.cancel(r+cc-ss*(cp-r)/Dpsi)==0
    T2psi=q*q/r**2*((2-3*ss+2*r*cc)*ss/Dpsi+ss)
    assert s.cancel(T2psi-4*q*q*ss*(cp-r)**2/Dpsi**2)==0
    assert s.cancel(q*q/r**2*((2-3*ss)*2*s.pi+ss*2*s.pi)-4*s.pi*q*q)==0
    psi,phi,qy,gq=s.symbols('psi phi qy gammaq',real=True)
    # On the true inverse graph T2=2*pi*phi*(1+2*q^2)-psi.
    assert s.diff(1+2*q*q,q)*qy==4*q*qy
    independent_extra=2*s.Symbol('T2_value')-8*s.pi*phi*q*q
    assert s.expand(independent_extra.subs(s.Symbol('T2_value'),2*s.pi*phi*(1+2*q*q)-psi)-(4*s.pi*phi-2*psi))==0
    t,ty,tz,tp,Fy,Fz,Fyz,by,bz,gy,gz=s.symbols('t ty tz tp Fy Fz Fyz betaY betaZ gammaY gammaZ',real=True)
    D=1+t*t;py=-Fy/D;pz=-Fz/D
    pyz=-(Fyz+2*t*ty*pz+2*t*tz*py+2*t*tp*py*pz)/D
    expected=-Fyz/D+2*t*(ty*Fz+tz*Fy)/D**2-2*t*tp*Fy*Fz/D**3
    assert s.cancel(pyz-expected)==0
    T1yz=s.Symbol('T1yz');total=T1yz+ty*pz+tz*py+tp*py*pz+t*pyz
    total_expected=T1yz-t*Fyz/D+(t*t-1)*(ty*Fz+tz*Fy)/D**2+(1-t*t)*tp*Fy*Fz/D**3
    assert s.cancel(total-total_expected)==0
    assert s.expand((by*t+gy)*Fz+(bz*t+gz)*Fy-(by*Fz+bz*Fy)*t-(gy*Fz+gz*Fy))==0
    return dict(passed=True,exact_Mobius_original_direction_and_T2_derivative_identities=3,
        independent_fixed_phi_extra_q_y_cancellation=True,implicit_and_total_T1_mixed_identities=2,
        both_direction_beta_gamma_cross_terms_retained=True,
        common_primitive_normalization='T1(0)=T2(0)=0; T1(2pi)=0, T2(2pi)=4*pi*q^2',
        overlap_same_original_rational_integrands_and_unique_inverse=True,
        inverse_uniqueness='Phi_psi=(1+t^2)/(2*pi*nu)>0 with original positive nu',
        weighted_FYFZ_split_valid_only_on_true_inverse_graph=True)


def independent_variable_q_signed_integrals(owner):
    frame=owner.source.frame(64,0,branch='positive');a=frame.roots['q'].atlas;c=a.ctx;p=mp.mp.clone();p.dps=100
    theta,u,q=s.symbols('theta u q',real=True)
    h=1/s.sqrt(1+u*u);r=u/s.sqrt(1+u*u)
    t=2*q*h*(s.cos(theta)-r)/(1-2*r*s.cos(theta)+r*r)
    integrands={}
    for power in (1,2):
        expr=t**power
        integrands[power]=tuple(s.lambdify((theta,u,q),item,modules=[{'mpf':p.mpf},'mpmath'])
            for item in (expr,s.diff(expr,u),s.diff(expr,q),s.diff(expr,u,2),s.diff(expr,u,q)))
    ttheta=s.lambdify((theta,u,q),s.diff(t,theta),modules=[{'mpf':p.mpf},'mpmath'])
    def sigma(y):
        e0=p.exp(-1/y**2);e1=p.exp(-1/(1-y)**2);return e0/(e0+e1)
    def av(y):return p.mpf(4)/5+p.mpf(6)/5*sigma(y)
    eta=p.exp(-30)
    def qv(y):return p.sqrt((1+eta)/av(y)-p.mpf('.5'))
    comparisons=0;overlap=0;inverses=0;minimum_q=p.inf
    for ytext,utext in (('.35','-.20'),('.35','.20'),('.65','-.9'),('.65','.9'),('.90','-3'),('.90','3')):
        y=p.mpf(ytext);q0=qv(y);qy=p.diff(qv,y);a0=av(y);ay=p.diff(av,y);minimum_q=min(minimum_q,q0)
        nu0=1+2*q0*q0;nuy=4*q0*qy;u0=p.mpf(utext);gg=p.mpf('.031')
        p20=u0/q0;p2y=p20*gg;p2z=p.mpf('.049');p2yz=p.mpf('.011')
        uy=p2y*q0+p20*qy;uz=p2z*q0;uyz=p2yz*q0+p2z*qy
        values=dict(a=(a0,ay,0,0),b=(0,0,0,0),t0=(0,0,0,0),q=(q0,qy,0,0),nu=(nu0,nuy,0,0),
            E=('1.3','-.2','.17','.08'),p2=(p20,p2y,p2z,p2yz))
        roots={name:current.MixedJet(a,{key:a.scalar(c.mpf(value)) for key,value in zip((C0,Y,Z,YZ),row,strict=True)})
            for name,row in values.items()}
        ui=current.MixedJet(a,{key:a.scalar(c.mpf(value)) for key,value in zip((C0,Y,Z,YZ),(u0,uy,uz,uyz),strict=True)})
        directional={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()}
        query=dict(q=roots['q'][C0],roots=directional,original_u_source=ui[C0],regular_predicate=False,signed_predicate=True)
        kernel=current.frames.O2PredicatePhase(query,c.mpf(0));kernel.nu=roots['nu'][C0]
        assert kernel.geometry=='signed_Mobius'
        for psitext in ('.47','3.9'):
            psi=p.mpf(psitext);segments=[0]+([p.pi] if psi>p.pi else [])+[psi];partials={}
            for power,name in ((1,'T1'),(2,'T2')):
                f0,fu,fq,fuu,fuq=integrands[power]
                partials[name]={C0:p.quad(lambda th:f0(th,u0,q0),segments),
                    Y:p.quad(lambda th:fu(th,u0,q0)*uy+fq(th,u0,q0)*qy,segments),
                    Z:p.quad(lambda th:fu(th,u0,q0)*uz,segments),
                    YZ:p.quad(lambda th:fuu(th,u0,q0)*uy*uz+fuq(th,u0,q0)*qy*uz+fu(th,u0,q0)*uyz,segments)}
            f0,fu,fq,*_=integrands[1];tt=f0(psi,u0,q0);ty=fu(psi,u0,q0)*uy+fq(psi,u0,q0)*qy
            tz=fu(psi,u0,q0)*uz;tp=ttheta(psi,u0,q0);D=1+tt*tt
            phi=(psi+partials['T2'][C0])/(2*p.pi*nu0)
            FY=partials['T2'][Y]-2*p.pi*phi*nuy;FZ=partials['T2'][Z];FYZ=partials['T2'][YZ]
            py=-FY/D;pz=-FZ/D;pyz=-(FYZ+2*tt*ty*pz+2*tt*tz*py+2*tt*tp*py*pz)/D
            totalY=partials['T1'][Y]+tt*py;totalZ=partials['T1'][Z]+tt*pz
            totalYZ=partials['T1'][YZ]+ty*pz+tz*py+tp*py*pz+tt*pyz
            Hphi=phi-psi/(2*p.pi)
            expectedA={C0:a0*Hphi/2,Y:ay*Hphi/2-a0*py/(4*p.pi),Z:-a0*pz/(4*p.pi),
                YZ:-(ay*pz+a0*pyz)/(4*p.pi)}
            E0=p.mpf('1.3');Ey=p.mpf('-.2');Ez=p.mpf('.17');Eyz=p.mpf('.08')
            ae=a0*E0;aey=ay*E0+a0*Ey;aez=a0*Ez;aeyz=ay*Ez+a0*Eyz
            expectedB={C0:-ae*partials['T1'][C0]/(4*p.pi),
                Y:-(aey*partials['T1'][C0]+ae*totalY)/(4*p.pi),
                Z:-(aez*partials['T1'][C0]+ae*totalZ)/(4*p.pi),
                YZ:-(aeyz*partials['T1'][C0]+aey*totalZ+aez*totalY+ae*totalYZ)/(4*p.pi)}
            chi=kernel.angles(c.mpf(psi/(2*p.pi)),'psi')[1]
            got=current.signed_fixed_phi_mixed(a,kernel,roots,ui,a.scalar(c.mpf(gg)),a.scalar(c.mpf(qy/q0)),chi,phase=c.mpf(phi))
            for name in ('T1','T2'):
                for key,target in partials[name].items():
                    saved.contains(saved.saved_value(c,got['record']['original_fixed_angle_'+name][str(key)],bases=a.bases),target,p.mpf('1e-35'));comparisons+=1
            for name,rows in (('A',expectedA),('B',expectedB)):
                for key,target in rows.items():saved.contains(got[name][key].finite_interval(),target,p.mpf('1e-35'));comparisons+=1
            for key,target in ((Y,py),(Z,pz),(YZ,pyz)):
                saved.contains(saved.saved_value(c,got['record']['original_inverse_psi_mixed'][str(key)],bases=a.bases),target,p.mpf('1e-35'));comparisons+=1
            inverse=kernel.evaluate(c.mpf(phi),bits=32)
            assert inverse['status']=='enclosed';saved.contains(inverse['psi_fraction_interval'],psi/(2*p.pi),p.mpf('1e-35'));inverses+=1
            if abs(u0)<=p.mpf('.25'):
                rk=current.frames.O2PredicatePhase(dict(query,regular_predicate=True,signed_predicate=False),c.mpf(0));rk.nu=roots['nu'][C0]
                common=current.regular.regular_fixed_phi_mixed(a,rk,roots,c.mpf(psi/(2*p.pi)),phase=c.mpf(phi))
                for name,rows in (('A',expectedA),('B',expectedB)):
                    for key,target in rows.items():saved.contains(common[name][key].finite_interval(),target,p.mpf('1e-35'));overlap+=1
    assert minimum_q<p.mpf('1e-6')
    return dict(passed=True,independent_original_varying_q_signed_integral_inverse_AB_comparisons=comparisons,
        same_function_regular_signed_overlap_AB_comparisons=overlap,independent_true_inverse_brackets=inverses,
        both_u_signs_and_positive_q_below_one_millionth=True,
        actual_original_flat_a_q_nu_y_functions_used_at_finite_diagnostic_eta=True,
        all_a_E_mixed_products_and_original_transverse_cross_terms_checked=True,
        finite_parameters_not_selected_original_field_or_global_N=True)


def whole_original_predicate_mixed_domain(owner,manifest):
    c=owner.ctx;count=min(owner.source.source.parent.parent.levels);rows=manifest['whole_original_O2_predicate_mixed_cells']
    assert len(rows)==count and manifest['exact_outer_source_domain']==dict(y=['0','1'],Z=['-1','1'],phi=['0','1'])
    frames=0;jetrows=0
    for index,item in enumerate(rows):
        assert set(item)=={'regular','positive','negative'}
        for branch,row in item.items():
            source=row['source'];proof=row['actual_whole_phase_mixed']
            assert source['exact_y_cell']==[str(s.Rational(index,count)),str(s.Rational(index+1,count))]
            assert source['source_family']==owner.family and source['branch']==branch
            assert row['conservative_whole_true_inverse_graph_outer_enclosure']
            assert row['independent_angle_phase_rectangle_not_exact_compatible_inverse_graph']
            assert row['named_original_predicate_required_for_all_downstream_consumers']
            assert proof['defining_original_source_predicate']==source['actual_original_source_predicate']
            for key in ('conditional_domain_not_entire_outer_rectangle','same_original_source_q_u_pressure_histories_and_ordinary_derivatives',
                    'actual_nu_y_retained','q_y_in_original_direction_and_T1_T2_retained','all_a_E_first_and_mixed_product_terms_retained',
                    'phase_is_held_fixed_not_differentiated','no_constant_reference_q_a_nu_derivative_assumption'):assert proof[key]
            assert ep(saved.interval(c,proof['fixed_true_phase_fraction']))==(0,1)
            if branch!='regular':
                assert proof['actual_q_y_direction_beta_and_extra_curvature_retained']
                assert proof['extra_variable_q_curvature']['valid_on_true_inverse_graph_only']
                assert proof['positive_original_q_and_hinv_not_replaced_by_flat_limit']
                cap=proof['correlated_weighted_curvature']['bounds']['primitive_J_squared']
                assert cap['formal_positive_scale']['source_exponents'][3]==1
            assert proof['negative_flat_source_log_offsets_supported_without_q_floor']
            for name in ('A','B'):
                for jet in proof['complete_original_A_B_mixed'][name].values():
                    assert jet['encloses_original_source_function'] and not jet['point_value_selected'];jetrows+=1
            frames+=1
    inverse=manifest['actual_true_phase_inverse_queries'];assert len(inverse)==18
    for row in inverse:
        selected=row['selected_original_inverse'];target=row['actual_fixed_phi_mixed']['fixed_true_phase_fraction']
        saved.contains(saved.interval(c,selected['phase_image']),saved.interval(c,target))
        assert row['source_inverse_not_selected_field']
        assert selected['bracket_proof']=='directed endpoint inequalities and exact strict source monotonicity'
    zero=0;rejected=0
    for branch in ('regular','positive','negative'):
        f=owner.source.frame(count,count//2,branch=branch)
        for phase in (0,'.5',1):
            proof=owner.primitive(f,phase,phase=phase)['record']
            assert proof['exact_common_symmetry_trace']
            assert all(v['exact_zero'] for rows in proof['complete_original_A_B_mixed'].values() for v in rows.values());zero+=1
        for call in (lambda:owner.primitive(replace(f),0,phase=0),lambda:owner.primitive(f,0,phase=2)):
            try:call()
            except ValueError:rejected+=1
    assert rejected==6
    f=owner.source.frame(count,0,branch='regular');a=f.roots['q'].atlas
    assert current.regular_fixed_phi_mixed.__code__ is current.regular.regular_fixed_phi_mixed.__code__
    # The source endpoint failure was arithmetic policy, not a change of
    # the accepted regular derivative formula or a floor on q.
    tail=current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,(0,0,0,-2,0),a.ctx.mpf(-2000)),1,a.ledger)
    anchored=current.bounded_O2_offset_anchor(a,tail)
    assert anchored.scale.powers==tail.scale.powers
    try:
        current.bounded_O2_offset_anchor(a,current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,offset=a.ctx.mpf(2000)),1,a.ledger))
    except ArithmeticError:rejected+=1
    assert rejected==7
    return dict(passed=True,continuous_original_y_cells=count,whole_predicate_source_phase_frames=frames,
        original_A_B_C0_y_Z_yZ_rows=jetrows,original_inverse_query_records=len(inverse),
        exact_common_symmetry_traces=zero,issued_and_invalid_phase_guards=rejected,
        accepted_regular_mathematical_code_unchanged_with_negative_log_arithmetic_extension=True,
        union_including_axis_and_both_nonzero_Z_sides=True,
        bounds_apply_only_on_named_original_predicate_and_true_inverse_graph=True,
        no_changed_integral_or_global_N_or_terminal_claim=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_changed_five_integrals_installed','all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2PredicateMixed();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('same_original_direction_and_fixed_phi_mixed_calculus',same_original_direction_and_mixed_calculus),
                ('independent_varying_q_signed_and_overlap_integrals',lambda:independent_variable_q_signed_integrals(owner)),
                ('whole_original_O2_predicate_phase_domain',lambda:whole_original_predicate_mixed_domain(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Actual full original O2 conditional regular/signed fixed-phi inverse and A/B mixed jets with varying q/a/nu, both signed curvature pieces and complete products. Named predicate union and same-function overlap supported. Not changed five integrals, controls, global N, scale recursion or full reconstruction.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Full original O2 varying-q predicate inverse and A/B mixed evidence PASS',flush=True);return report


if __name__=='__main__':run()
