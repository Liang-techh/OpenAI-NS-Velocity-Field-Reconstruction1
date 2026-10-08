"""Independent original variable-q integral and fixed-phi mixed evidence."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_regular_mixed_phase as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

C0,Y,Z,YZ=current.C0,current.Y,current.Z,current.YZ;ep=current.ep


def fixed_phi_calculus():
    q,qy,P,Py,Pz,Pyz,phi=s.symbols('q qy P Py Pz Pyz phi',real=True)
    nu=1+2*q*q;nuy=4*q*qy
    FY=q*q*Py+2*q*qy*P-2*s.pi*phi*nuy
    assert s.expand(FY-(q*q*Py+nuy/2*(P-4*s.pi*phi)))==0
    assert s.expand(q*q*Pyz+2*q*qy*Pz-(q*q*Pyz+nuy/2*Pz))==0
    t,ty,tz,tp,Fy,Fz,Fyz=s.symbols('t ty tz tp Fy Fz Fyz',real=True)
    D=1+t*t;py=-Fy/D;pz=-Fz/D
    pyz=-(Fyz+2*t*ty*pz+2*t*tz*py+2*t*tp*py*pz)/D
    expanded=-Fyz/D+2*t*(ty*Fz+tz*Fy)/D**2-2*t*tp*Fy*Fz/D**3
    assert s.cancel(pyz-expanded)==0
    T1yz=s.Symbol('T1yz');total=T1yz+ty*pz+tz*py+tp*py*pz+t*pyz
    expected=T1yz-t*Fyz/D+(t*t-1)*(ty*Fz+tz*Fy)/D**2+(1-t*t)*tp*Fy*Fz/D**3
    assert s.cancel(total-expected)==0
    return dict(passed=True,exact_source_q_nu_cancellations=2,
        independent_implicit_and_total_primitive_chain_rule_identities=2,
        inverse_uniqueness='dPhi/dpsi=(1+t^2)/(2*pi*nu)>0; original nu>0',
        phase_normalization_is_not_NS_viscosity=True)


def independent_variable_q_integrals(owner):
    count=min(owner.source.parent.parent.levels);frame=owner.source.source_frame(count,0,Z_lower=0,Z_upper=0)
    a=frame.roots['q'].atlas;c=a.ctx;p=mp.mp.clone();p.dps=100
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
    eta=p.exp(-12)
    def qv(y):return p.sqrt((1+eta)/av(y)-p.mpf('.5'))
    comparisons=0;minimum_q=p.inf
    for ytext,utext in (('.35','-.12'),('.65','.12'),('.90','-.12')):
        y=p.mpf(ytext);q0=qv(y);qy=p.diff(qv,y);a0=av(y);ay=p.diff(av,y);minimum_q=min(minimum_q,q0)
        nu0=1+2*q0*q0;nuy=4*q0*qy;u0=p.mpf(utext)
        p20=u0/q0;p2y=p.mpf('.031');p2z=p.mpf('.049');p2yz=p.mpf('.011')
        uy=p2y*q0+p20*qy;uz=p2z*q0;uyz=p2yz*q0+p2z*qy
        values=dict(a=(a0,ay,0,0),b=(0,0,0,0),t0=(0,0,0,0),q=(q0,qy,0,0),nu=(nu0,nuy,0,0),
            E=('1.3','-.2','.17','.08'),p2=(p20,p2y,p2z,p2yz))
        roots={name:current.MixedJet(a,{key:a.scalar(c.mpf(value)) for key,value in zip((C0,Y,Z,YZ),row,strict=True)})
            for name,row in values.items()}
        directional={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()}
        kernel=current.source.positive.PositiveLogQPhase(dict(q=roots['q'][C0],roots=directional,
            original_u_source=roots['p2'][C0]*roots['q'][C0]),c.mpf(0))
        assert kernel.geometry=='small_r_series'
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
            got=current.regular_fixed_phi_mixed(a,kernel,roots,c.mpf(psi/(2*p.pi)),phase=c.mpf(phi))
            for name in ('T1','T2'):
                for key,target in partials[name].items():
                    saved.contains(saved.saved_value(c,got['record']['original_fixed_angle_'+name][str(key)],bases=a.bases),target,p.mpf('1e-36'));comparisons+=1
            for name,rows in (('A',expectedA),('B',expectedB)):
                for key,target in rows.items():
                    saved.contains(got[name][key].finite_interval(),target,p.mpf('1e-36'));comparisons+=1
            for key,target in ((Y,py),(Z,pz),(YZ,pyz)):
                saved.contains(saved.saved_value(c,got['record']['original_inverse_psi_mixed'][str(key)],bases=a.bases),target,p.mpf('1e-36'));comparisons+=1
    assert minimum_q<p.mpf('.01')
    return dict(passed=True,independent_original_variable_q_T1_T2_inverse_A_B_mixed_comparisons=comparisons,
        both_regular_u_signs_and_q_below_one_half_and_one_hundredth=True,
        actual_original_flat_profile_a_q_y_functions_used_at_finite_diagnostic_eta=True,
        full_fixed_phi_nu_y_and_a_E_product_terms_checked=True,
        finite_parameters_are_diagnostics_only_not_original_selected_data=True)


def original_axis_domain(owner,manifest):
    count=min(owner.source.parent.parent.levels);rows=manifest['whole_original_axis_source_cells']
    assert manifest['exact_source_domain']==dict(y=['0','1'],Z=['0','0'],phi=['0','1'])
    assert len(rows)==count;jetrows=0
    for index,item in enumerate(rows):
        frame=item['source'];proof=item['actual_whole_phase_mixed']
        assert item['conservative_whole_phase_inverse_graph_outer_enclosure']
        assert item['independent_angle_phase_rectangle_not_exact_compatible_inverse_graph']
        assert item['accepted_original_q_nu_identity_receipt']==current.source.RECEIPT
        assert frame['source_family']==owner.family
        assert frame['exact_y_cell']==[str(s.Rational(index,count)),str(s.Rational(index+1,count))]
        assert frame['exact_Z_range']==['0','0']
        assert frame['actual_root_jets']['p2'][str(C0)]['exact_zero']
        assert not frame['actual_root_jets']['q'][str(C0)]['exact_zero']
        assert not frame['actual_root_jets']['q'][str(Y)]['exact_zero']
        assert not frame['actual_root_jets']['nu'][str(Y)]['exact_zero']
        assert proof['actual_nu_y_retained'] and proof['q_y_in_original_direction_and_T1_T2_retained']
        assert proof['phase_is_held_fixed_not_differentiated'] and proof['all_a_E_first_and_mixed_product_terms_retained']
        assert proof['no_constant_reference_q_a_nu_derivative_assumption'] and proof['no_division_by_r_or_p2']
        assert ep(saved.interval(owner.ctx,proof['fixed_true_phase_fraction']))==(0,1)
        assert ep(saved.interval(owner.ctx,proof['source_angle_fraction']))==(0,1)
        for name in ('A','B'):
            for key,jet in proof['complete_original_A_B_mixed'][name].items():
                assert jet['encloses_original_source_function'] and not jet['point_value_selected'];jetrows+=1
                assert jet['formal_positive_scale']['radius_power']==0
    inverse=manifest['actual_true_phase_inverse_queries'];assert len(inverse)==6
    for item in inverse:
        selected=item['selected_original_inverse'];target=item['actual_fixed_phi_mixed']['fixed_true_phase_fraction']
        assert item['accepted_original_q_nu_identity_receipt']==current.source.RECEIPT
        assert selected['chart']=='psi' and item['source_inverse_not_selected_field']
        saved.contains(saved.interval(owner.ctx,selected['phase_image']),saved.interval(owner.ctx,target))
        assert selected['bracket_proof']=='directed endpoint inequalities and exact strict source monotonicity'
    # Fresh symmetry and unsupported-geometry guards only.
    frame=owner.source.source_frame(count,count//2,Z_lower=0,Z_upper=0);a=frame.roots['q'].atlas
    proof=owner.evaluate(frame,'.5')['actual_fixed_phi_mixed']
    assert proof['exact_common_symmetry_trace'] and all(v['exact_zero'] for row in proof['complete_original_A_B_mixed'].values() for v in row.values())
    rejected=0
    for call in (lambda:current.regular_fixed_phi_mixed(a,owner.kernel(frame),frame.roots,0,phase=2),
        lambda:owner.kernel(owner.source.source_frame(count,count//2,Z_lower=-1,Z_upper=1))):
        try:call()
        except (ValueError,ArithmeticError):rejected+=1
    assert rejected==2
    return dict(passed=True,continuous_original_axis_source_and_whole_true_phase_cells=count,
        ordinary_y_Z_yZ_AB_jet_rows_checked=jetrows,source_inverse_bracket_records=6,
        exact_common_halfperiod_trace_and_unsupported_guards=3,
        full_original_pressure_P0_and_positive_q_source_retained=True,
        no_neighboring_Z_signed_or_changed_five_integral_claim=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('whole_neighboring_Z_or_signed_domain_coverage','actual_changed_five_integrals_installed',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed','current_whole_N_selected',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2RegularMixed();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('fixed_phi_original_source_calculus',fixed_phi_calculus),
            ('independent_variable_q_defining_integrals',lambda:independent_variable_q_integrals(owner)),
            ('whole_original_O2_axis_phase_domain',lambda:original_axis_domain(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Actual varying-q original O2 regular fixed-phi inverse and full A/B C0/y/Z/yZ, entire y[0,1] at Z0, source-backed ordinary transverse jets. Not nonzero-Z/signed coverage, changed five integrals or full reconstruction.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Actual O2 varying-q regular mixed inverse and full A/B source evidence PASS',flush=True)
    return report


if __name__=='__main__':run()
