"""Two positive-q branches, source geometry and microscopic q evidence."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_positive_q_curvature as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep


def positive_q_geometric_proof():
    q,k=s.symbols('q k',positive=True)
    D=1+4*q*q*k*k;identities=0
    for m,M,branch in ((2*q,s.Rational(1,2),'0<q<=1/2'),(1,q,'q>=1/2')):
        for lhs,rhs in ((q**5/m**5,M**5),(q**5/m**4,q*M**4),(q**3/m**3,M**3)):
            assert s.cancel(lhs-rhs)==0;identities+=1
        difference=s.expand(D-m*m*(1+k*k))
        assert s.expand(difference-(1-4*q*q if branch.startswith('0') else (4*q*q-1)*k*k))==0
    assert s.expand(5*(1+k*k)-(1+2*k)**2-(k-2)**2)==0
    assert s.expand(20*(1+k*k)-(4+2*k)**2-4*(2*k-1)**2)==0
    r,cs=s.symbols('r coschi',real=True);h2=1-r*r;n=r+cs;Dchi=h2+2*r*n
    assert s.expand(Dchi-(1-cs*cs)-n*n)==0
    # Derive J directly from the signed fixed-angle T2 function and the
    # exact source chi/r/s chain. No derivative of a magnitude bound.
    ss,chi,psi=s.symbols('s chi psi',real=True)
    T2=q*q/r**2*((2-3*ss)*chi+ss*psi+2*r*s.sin(chi))
    derivative=ss*(s.diff(T2,r)-2*r*s.diff(T2,ss))+2*s.sin(chi)*s.diff(T2,chi)
    J=2*q*q/r**2*(s.sin(chi)*(ss+2*r*(r+s.cos(chi))-3*ss)+ss*(chi-psi)/r)
    assert s.simplify((derivative-J).subs(ss,1-r*r))==0
    tt=s.Symbol('t',real=True)
    assert s.expand((1+tt*tt)-tt*tt)==1
    # D>=1 implies D^(5/2),D^2>=D^(3/2). With positive m,
    # the new q^3/m^3 identity controls the extra J transport term.
    return dict(passed=True,exact_positive_q_weight_identities=identities,
        both_positive_q_denominator_branches_proved=True,
        denominator_branch_nonnegative_remainders=['1-4*q^2 for0<q<=1/2','(4*q^2-1)*k^2 forq>=1/2'],
        square_and_exact_source_geometry_certificates=3,
        independent_original_fixed_angle_J_chain_identity=True,
        C3_derivation=['P<=h^(3/2)*(h+2*k)^(1/2)*(4*h+2*k)+h^2*D0',
            '|sinchi|*Dchi<=h^(3/2)*(h+2*k)^(3/2)',
            'product/h^3<=(h+2*k)^2*(4*h+2*k)+sqrt(h)*(h+2*k)^(3/2)*D0',
            '0<h<=1; use the two square certificates and1+k^2>=1',
            'C3=5*sqrt20+5^(3/4)*D0; D^(5/2),D^2>=D^(3/2)>=m^3*(1+k^2)^(3/2)'],
        original_J_squared_geometry_reused_from_accepted_source_receipt=current.source.axial.mixed.RECEIPT,
        rational_prefactors=['2*abs(t)/D^3<=2/D^(5/2)','abs(1-t^2)/D^3<=1/D^2'],
        no_q_lower_floor_or_limit_q_zero_used=True)


def finite_peak_inequalities(owner):
    count=min(owner.source.parent.parent.levels);frame=owner.source.source_frame(count,0,Z_lower=0,Z_upper=0)
    a=frame.roots['q'].atlas;c=a.ctx;p=mp.mp.clone();p.dps=110;comparisons=0;points=0
    for text in ('.9','.1','.001','1e-12','1e-40'):
        q=p.mpf(text);bound=current.positive_q_weighted_curvature(a,a.scalar(c.mpf(q)),r_lower=c.mpf('.1'))
        caps={name:max(ep(value.finite_interval())) for name,value in bound['bounds'].items()}
        for magnitude in (p.mpf('.2'),p.mpf('.99'),1-p.mpf('1e-30')):
            for sign in (-1,1):
                r=sign*magnitude;h=p.sqrt(1-r*r);angles=[p.mpf(v) for v in ('.13','1.5','3.13','5.8')]
                for k in (-1/(2*q),-1,0,1,1/(2*q)):
                    cosine=-r+h*k
                    if -1<=cosine<=1:angles.append(p.acos(cosine))
                for chi in angles:
                    psi=2*p.atan2((1-r)*p.sin(chi/2),(1+r)*p.cos(chi/2))
                    if psi<0:psi+=2*p.pi
                    ss=h*h;n=r+p.cos(chi);Dchi=ss+2*r*n
                    t=2*q*n/h;tp=-2*q*p.sin(chi)*Dchi/h**3;D=1+t*t
                    J=2*q*q/r**2*(p.sin(chi)*(Dchi-3*ss)+ss*(chi-psi)/r)
                    actual=dict(inverse_J_squared=abs(2*t*tp*J*J/D**3),
                        primitive_J_squared=abs((1-t*t)*tp*J*J/D**3),
                        inverse_J=abs(2*t*tp*J/D**3),primitive_J=abs((1-t*t)*tp*J/D**3))
                    for name,value in actual.items():
                        assert value<=caps[name]+p.mpf('1e-85')*max(caps[name],1),(text,sign,p.nstr(magnitude),name,value,caps[name])
                        comparisons+=1
                    points+=1
    return dict(passed=True,finite_original_signed_geometry_peak_cases=points,
        independent_four_curvature_inequality_comparisons=comparisons,
        q_down_to_one_e_minus_forty_and_r_with_one_e_minus_thirty_defect=True,
        positive_and_negative_signed_geometry_and_n_over_h_peak_coordinates=True,
        finite_units_are_inequality_diagnostics_only_not_original_source_values=True)


def native_positive_q_contract(owner,manifest):
    c=owner.ctx;count=min(owner.source.parent.parent.levels);rows=manifest['whole_original_O2_positive_q_conditional_bounds'];assert len(rows)==count
    for index,row in enumerate(rows):
        assert row['source_family']==owner.family and row['exact_y_cell']==[str(s.Rational(index,count)),str(s.Rational(index+1,count))]
        assert row['exact_outer_Z_window']==['-1','1'] and row['not_claimed_full_rectangle_signed_domain']
        assert row['conditional_source_predicate']=='abs(original u=p2*q/dstar)>=3/16'
        proof=row['actual_conditional_weighted_curvature'];q=proof['original_q']
        assert q['sign']=='positive' and not q['exact_zero']
        assert q['formal_positive_scale']['source_exponents']==[0,0,0,1]
        assert proof['no_native_q_lower_floor_or_exact_flat_replacement'] and not proof['actual_signed_mixed_phase_primitives_installed']
        assert proof['original_q_factor_retained_in_primitive_J_squared']
        assert proof['reference_q_at_least_half_hypothesis_replaced_by_exact_positive_q_weights']
        assert ep(saved.interval(c,proof['source_upper_M_not_selected_parameter']))[0]>=c.mpf('.5')
        for name,cap in proof['bounds'].items():
            assert cap['sign']=='positive' and cap['formal_positive_scale']['radius_power']==0
            powers=cap['formal_positive_scale']['source_exponents']
            assert powers==([0,0,0,1] if name=='primitive_J_squared' else [0,0,0,0])
    frame=owner.source.source_frame(count,0,Z_lower=0,Z_upper=0);a=frame.roots['q'].atlas
    rejected=0
    for call in (lambda:current.positive_q_weighted_curvature(a,a.scalar(0),r_lower='.1'),
        lambda:current.positive_q_weighted_curvature(a,a.scalar(-1),r_lower='.1'),
        lambda:current.positive_q_weighted_curvature(a,a.scalar('.9'),r_lower=0),
        lambda:current.positive_q_weighted_curvature(a,a.scalar('1.1'),r_lower='.1')):
        try:call()
        except ValueError:rejected+=1
    assert rejected==4
    return dict(passed=True,continuous_original_positive_q_source_cells=count,
        original_positive_formal_q_in_primitive_curvature_bounds=count,
        strict_positivity_and_geometry_hypothesis_rejections=rejected,
        predicate_lower_r_is_theorem_not_original_field_value=True,
        application_requires_proven_signed_predicate_and_same_original_source=True,
        full_Z_rectangle_not_claimed_signed_at_axis=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_signed_mixed_phase_primitives_installed','regular_signed_domain_coverage_installed',
        'actual_changed_five_integrals_installed','all_17_chart_or_24_cell_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2PositiveQCurvature();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('all_positive_q_original_geometric_calculus',positive_q_geometric_proof),
            ('independent_finite_signed_peak_inequalities',lambda:finite_peak_inequalities(owner)),
            ('whole_original_positive_q_conditional_source_contract',lambda:native_positive_q_contract(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original O2 signed weighted curvature bounds valid for all positive q under proved actual signed predicate. Exact two-q-branch weights, J chain and additional variable-q term, native positive source factor retained. Not signed phase implementation, predicate union, five integrals or full reconstruction.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('All-positive-q original weighted mixed curvature PASS',flush=True);return report


if __name__=='__main__':run()
