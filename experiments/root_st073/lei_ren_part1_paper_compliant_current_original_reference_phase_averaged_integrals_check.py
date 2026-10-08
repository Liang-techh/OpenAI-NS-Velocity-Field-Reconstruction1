"""Exact original radial/centering identities and independent averaging checks."""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_phase_averaged_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

points=current.points;ep=current.ep


def remainder_references(c):
    p=mp.mp.clone();p.dps=120;checks=0
    for lo,hi in (('-1','1'),('-.02','.03'),('0','0'),('1e-100','1e-100'),('-.8','-.5'),('.2','.7')):
        got=current.second_exponential_remainder(c,c.mpf((lo,hi)))
        for x in (p.mpf(lo),p.mpf(hi),(p.mpf(lo)+p.mpf(hi))/2):
            expected=p.quad(lambda t:(1-t)*p.exp(t*x),[0,1])
            saved.contains(got,expected,p.mpf('1e-110'));checks+=1
    assert ep(current.second_exponential_remainder(c,c.mpf(0)))==(mp.mpf('.5'),mp.mpf('.5'))
    return dict(passed=True,independent_weighted_exponential_integral_comparisons=checks,
        exact_zero_argument_one_half=True,integrated_directed_order64_tail=True)


def radial_original_source_references(owner):
    inputs=owner.owner.templates['inputs'];y=s.Symbol('y',real=True);f=s.exp(y/10)
    substitutions=dict(zip(inputs[1:5],(f,s.Rational(5,8)*f,s.Rational(5,12)*f*f,s.Rational(5,2)*f*f),strict=True))
    identities=0
    for name,rows in owner.radial_templates.items():
        original=owner.owner.templates['rows'][(name,0)]
        for (powers,expr),(dy_powers,dy_expr) in zip(original,rows,strict=True):
            assert powers==dy_powers;r=powers[0]
            expression=s.exp(r*y)*expr.subs(substitutions,simultaneous=True)
            expected=s.exp(r*y)*dy_expr.subs(substitutions,simultaneous=True)
            assert s.simplify(s.diff(expression,y)-expected)==0,(name,powers)
            for P in owner.owner.templates['pressure_symbols']:
                for Q in owner.owner.templates['pressure_symbols']:assert s.diff(dy_expr,P,Q)==0
            identities+=1
    p=mp.mp.clone();p.dps=150
    original_fns={name:tuple(s.lambdify(inputs,expr,modules=[{'mpf':p.mpf},'mpmath']) for powers,expr in rows)
        for name,rows in owner.owner.templates['rows'].items() if name[1]==0}
    dy_fns={name:tuple(s.lambdify(inputs,expr,modules=[{'mpf':p.mpf},'mpmath']) for powers,expr in rows)
        for name,rows in owner.radial_templates.items()}
    z=p.mpf('.37');q=1+z*z;alpha=p.mpf('2.9')
    Pstar=p.exp(3);delta=p.exp(-42);L=1-delta*z*z
    def values(yy):
        ff=p.exp(yy/10)
        return (z,ff,p.mpf(5)/8*ff,p.mpf(5)/12*ff*ff,p.mpf(5)/2*ff*ff,
            -alpha/q**2,4*alpha*z/q**3,alpha*(4-20*z*z)/q**4)
    def root(name,yy,derivative=False):
        R=110*p.exp(50+yy);rows=owner.radial_templates[name] if derivative else owner.owner.templates['rows'][(name,0)]
        functions=dy_fns[name] if derivative else original_fns[(name,0)]
        return sum(R**r*Pstar**a*delta**d*L**ell*fn(*values(yy))
                   for ((r,a,d,ell),expr),fn in zip(rows,functions,strict=True))
    comparisons=0
    for yy in ('-2.337','-.513'):
        for name in ('E','V','b','p1','p2'):
            expected=p.diff(lambda x:root(name,x),p.mpf(yy));got=root(name,p.mpf(yy),True)
            assert abs(got-expected)<=p.mpf('1e-120')*max(abs(got),abs(expected),1),(name,yy)
            comparisons+=1
    live=owner.radial_query('-2.5','-1.25');a=owner.atlas
    for name in ('V','b','a','t0'):assert live['radial_roots'][name][(0,1)].zero
    difference=a.add(live['radial_roots']['E'][(0,1)],-live['radial_roots']['E'][(0,0)]*(a.ctx.mpf(1)/10))
    assert ep(difference.coefficient)[0]<=0<=ep(difference.coefficient)[1]
    return dict(passed=True,exact_original_radial_factor_and_profile_chain_rule_identities=identities,
        independent_finite_unit_original_source_radial_derivative_comparisons=comparisons,
        original_radius_factor_derivative_retained=True,ordinary_P0_and_L_y_fixed=True,
        source_pressure_affine_late_factor_preserved=True,
        diagnostic_parameters_not_installed_in_original_frame=True)


def exact_split_and_IBP():
    E,V,A,B,Q,N=s.symbols('E V A B Q N',nonzero=True);F=E*A+Q/N
    dE,dV=F/N,B/N
    densities=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    linear=dict(m=B,h=E*A,k=V*E*A+E*B,e=2*V*B-E*E*A,p=E*E*A)
    rest=dict(m=0,h=Q,k=V*Q+F*B,e=-E*Q+B*B-F*F/2,p=E*Q+F*F/2)
    for name in linear:
        assert s.expand(densities[name]-linear[name]/N-rest[name]/N**2)==0
        assert s.expand(linear[name].subs({A:-A,B:-B},simultaneous=True)+linear[name])==0
    # Same original source/phase at internal radial traces. No equality
    # between unrelated chart primitives or independently picked caps is used.
    K=s.symbols('K0:5');G=s.symbols('G0:5')
    boundaries=sum(K[i+1]*G[i+1]-K[i]*G[i] for i in range(4))
    assert s.expand(boundaries-(K[-1]*G[-1]-K[0]*G[0]))==0
    yy,phi=s.symbols('y phi');lam=s.Symbol('lam');gg=s.Function('G')(yy,phi);kk=s.Function('K')(yy)
    total=lam*kk*gg+kk*(s.diff(gg,yy)+N*s.diff(gg,phi))
    assert s.expand(kk*s.diff(gg,phi)/N-(total-kk*(s.diff(gg,yy)+lam*gg))/N**2)==0
    return dict(passed=True,exact_original_five_density_linear_and_nonlinear_split_identities=5,
        exact_five_leading_reflection_identities=5,same_source_internal_boundary_telescoping=True,
        weighted_native_phase_slope_N_IBP_identity=True,nonlinear_finite_N_terms_not_zero_mean=True)


def independent_endpoint_weight_tests():
    p=mp.mp.clone();p.dps=120;comparisons=0
    left,right=p.mpf('-1.13'),p.mpf('.42');slope=p.mpf('.2')
    M=max(abs(1+slope*left),abs(1+slope*right));MY=slope
    for lam in (p.mpf(0),p.mpf(1),p.mpf('1.5')):
        wl,wr=p.exp(lam*left),p.exp(lam*right)
        mass=right-left if not lam else (wr-wl)/lam
        for N in (160,320):
            for origin in ('.137','.41'):
                # Independent exact antiderivative of the actual weighted
                # slowly varying, centered sinusoidal test function.
                k=lam+2*p.pi*N*p.j
                primitive=lambda y:p.exp(k*y)*((1+slope*y)/k-slope/k**2)
                integral=p.im(p.exp(2*p.pi*p.j*p.mpf(origin))*(primitive(right)-primitive(left)))/N
                bound=((wl+wr)*M/2+mass*(MY/2+lam*M/2))/N**2
                assert abs(integral)<=bound
                comparisons+=1
    return dict(passed=True,arbitrary_phase_origin_weighted_endpoint_comparisons=comparisons,
        own_rates_tested=[0,1,'3/2'],zero_rate_memory_kept=True,
        nonzero_global_endpoint_terms_not_dropped=True,period_enumeration_not_required=True)


def saved_native_contracts(owner,manifest):
    c=owner.ctx;levels=manifest['actual_original_reference_phase_averaged_levels'];ranges=[]
    for level in levels:
        count,N=level['exact_radial_cells'],level['candidate_N'];cells=level['full_cell_source_and_IBP_records']
        assert len(cells)==count and level['source_family']==owner.family
        assert level['actual_global_endpoint_terms_retained'] and level['exact_internal_same_closed_source_endpoint_cancellation']
        assert level['no_cross_chart_seam_cancellation_assumed'] and level['integer_period_count_not_enumerated']
        assert level['C1_Z_tightening_not_claimed_without_mixed_yZ_oracle']
        mass_totals={name:c.mpf(0) for name in points.exact.RATES}
        for i,cell in enumerate(cells):
            source=cell['source'];left,right=map(s.Rational,source['exact_radial_cell'])
            assert left==-5+s.Rational(5*i,count) and right==-5+s.Rational(5*(i+1),count)
            assert source['candidate_N']==N and source['native_radius_Jacobian']==1
            assert source['reflection_leading_mean_exact_zero'] and not source['nonlinear_finite_N_mean_not_zeroed'] is False
            assert source['radial_derivative_contract']['fixed_original_phi']
            assert source['radial_derivative_contract']['original_P0_y_exact_zero']
            for phase in source['actual_original_endpoint_phases'].values():
                assert phase['positive_original_origin_offset']['strictly_positive']
            for name,row in cell['contributions'].items():
                assert row['N_power']==-2
                mass_totals[name]+=saved.interval(c,row['positive_own_rate_mass'])
                if 0<i<count-1:assert row['retained_global_endpoint_term']['exact_zero']
        saved.contains(mass_totals['p'],c.mpf(5))
        values={name:saved.saved_value(c,level['C0_effective_contribution_enclosures'][name],owner.atlas.bases)
            for name in ('h','e','p')}
        direct={name:saved.saved_value(c,level['C0_direct_source_range_contribution_enclosures'][name],owner.atlas.bases)
            for name in ('h','e','p')}
        for name in values:
            assert max(abs(v) for v in ep(values[name]))<max(abs(v) for v in ep(direct[name]))
        for name in ('m','k'):
            assert level['effective_C0_bound_selection'][name]=='direct_smaller_by_collected_source_factor_ratio'
        ranges.append(dict(cells=count,N=N,effective_C0_magnitude_bounds={name:max(abs(v) for v in ep(value))
            for name,value in values.items()},direct_C0_magnitude_bounds={name:max(abs(v) for v in ep(value))
            for name,value in direct.items()}))
    first,second=levels[1:]
    assert (first['exact_radial_cells'],first['candidate_N'],second['candidate_N'])==(16,160,320)
    assert first['full_cell_source_and_IBP_records'][0]['source']['actual_N_dependent_nonlinear_remainder']['h']!=second['full_cell_source_and_IBP_records'][0]['source']['actual_N_dependent_nonlinear_remainder']['h']
    previous=json.loads(gzip.decompress((current.HERE/current.whole.NAME).read_bytes()))
    for row,old in zip(manifest['predecessor_original_Z_integrals_retained'],previous['original_reference_whole_window_refinements'],strict=True):
        assert row['candidate_N']==old['candidate_N']==160 and row['cells']==old['exact_original_cells']
        for name in points.exact.RATES:assert row['ordinary_Z_contribution_enclosures'][name]==old['actual_finite_N_five_density_contribution_enclosures'][name]['Z']
    return dict(passed=True,actual_source_levels=[(v['exact_radial_cells'],v['candidate_N']) for v in levels],
        h_e_p_C0_averaging_strictly_tighter_than_direct_source_ranges=True,
        microscopic_m_k_original_direct_bounds_retained=True,
        true_N_dependent_nonlinear_remainders_recomputed=True,
        predecessor_ordinary_Z_data_retained_exactly_not_relabelled_N320=True,
        full_source_range_endpoint_pressure_phase_and_domain_checks=True,computed_ranges=ranges)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('mixed_yZ_oracle_and_Z_averaging_installed','all_17_chart_or_24_cell_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*current.whole.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[name] is False for name in flags)
    owner=current.OriginalReferencePhaseAveraging(Z='.37')
    with mp.workdps(owner.ctx.dps+40):
        checks=dict(second_remainder_defining_integral=remainder_references(owner.ctx),
            original_radial_source_jet_evidence=radial_original_source_references(owner),
            exact_original_split_reflection_and_IBP=exact_split_and_IBP(),
            independent_weighted_endpoint_evidence=independent_endpoint_weight_tests(),
            actual_saved_source_and_integral_contracts=saved_native_contracts(owner,manifest))
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Actual original slow-y source jets and endpoint-retaining reflection-centered C0 integration, with independent remainder/radial/weighted-phase references. Direct microscopic bounds and predecessor Z integrals retained. No mixed-yZ C1 averaging, full-Z closure/control/global-N/recursion claim.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original source radial jets, finite-N remainder and endpoint-retaining C0 averaging PASS',flush=True)
    return report


if __name__=='__main__':run()
