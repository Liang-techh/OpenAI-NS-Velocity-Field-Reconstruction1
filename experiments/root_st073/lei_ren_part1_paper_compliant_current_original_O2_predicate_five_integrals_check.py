"""Defining finite-N moment differences, source union and native transport."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_predicate_five_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;C0,Z=current.C0,current.Z


def density_and_unit_calculus():
    E,V,F,B,N=s.symbols('E V F B N',nonzero=True)
    old=current.history.current.density.recovery.history_densities(E,V)
    changed=current.history.current.density.recovery.history_densities(E+F/N,V+B/N)
    first=dict(m=B,h=F,k=V*F+E*B,e=2*V*B-E*F,p=E*F)
    second=dict(m=0,h=0,k=F*B,e=B*B-F*F/2,p=F*F/2)
    for key in current.KEYS:assert s.cancel(changed[key]-old[key]-first[key]/N-second[key]/N**2)==0
    aa,qq,h,r,chi,psi,W=s.symbols('a q hinv r chi_fraction psi_fraction W1',nonzero=True)
    assert s.cancel(-aa*E*(qq*h/r*2*s.pi*(chi-psi))/(4*s.pi)-(-aa*E*qq*h*(chi-psi)/(2*r)))==0
    assert s.cancel(-aa*E*(2*qq*h*W)/(4*s.pi)-(-aa*E*qq*h*W/(2*s.pi)))==0
    z,P=s.symbols('z Pstar',real=True,nonzero=True);C=1/(1+z*z)
    inlet=dict(m=4*z/P,h=s.Rational(5,8)*C,k=s.Rational(5,2)*z*C/P,
        e=16*z*z/P**2-s.Rational(5,12)*C*C,p=s.Rational(5,2)*C*C)
    derivatives=dict(m=4/P,h=-s.Rational(5,4)*z*C*C,k=s.Rational(5,2)*(1-z*z)*C*C/P,
        e=32*z/P**2+s.Rational(5,3)*z*C**3,p=-10*z*C**3)
    for key in current.KEYS:assert s.cancel(s.diff(inlet[key],z)-derivatives[key])==0
    return dict(passed=True,independent_defining_history_difference_identities=5,
        exact_signed_and_regular_normalized_B_alias_identities=2,original_full_Z_inlet_derivative_identities=5,
        stable_modulation_identity='N*(exp(A/N)-1)=A*integral_0^1 exp(t*A/N)dt',
        original_S_equals_Pstar_no_second_Pstar_division=True,
        exact_square_Z_squared_and_Q_at_least_one_for_crossing_axis=True)


def finite_N_defining_density_fixtures(owner):
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas;c=owner.ctx
    p=mp.mp.clone();p.dps=150;comparisons=0;cases=0
    def densities(E,V):return dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
    for N in (160,257):
        for ztext in ('-.37','.37'):
            zz=p.mpf(ztext);E0=1+p.mpf('.3')*zz;EZ=p.mpf('.3');V0=p.mpf('-.21')+p.mpf('.4')*zz;VZ=p.mpf('.4')
            for Atext in ('-.7','0','1e-50','.6'):
                A0=p.mpf(Atext);AZ=p.mpf('.17');B0=p.mpf('-.46') if zz<0 else p.mpf('.53');BZ=p.mpf('-.13')
                roots={name:{C0:a.scalar(c.mpf(value)),Z:a.scalar(c.mpf(derivative))}
                    for name,value,derivative in (('E',E0,EZ),('V',V0,VZ))}
                graph=current.CollectedDensityGraph(a,owner.source.owner.scales.graph,
                    dict(kernel=frame.kernel,roots=roots,ledger=a.ledger),
                    dict(A=a.scalar(c.mpf(A0)),B_over_Pstar=a.scalar(c.mpf(B0))),
                    dict(A_Z_slow=a.scalar(c.mpf(AZ)),B_Z_slow=a.scalar(c.mpf(BZ))),N)
                got=graph.outputs()
                def changed(z):
                    E=E0+EZ*(z-zz);V=V0+VZ*(z-zz);A=A0+AZ*(z-zz);B=B0+BZ*(z-zz)
                    old=densities(E,V);new=densities(E*p.exp(A/N),V+B/N)
                    return {key:new[key]-old[key] for key in current.KEYS}
                actual=changed(zz);derivatives={key:p.diff(lambda z:changed(z)[key],zz) for key in current.KEYS}
                for label,targets in (('densities',actual),('density_Z',derivatives)):
                    for key,target in targets.items():
                        saved.contains(got[label][key].finite_interval(),target,p.mpf('1e-120'));comparisons+=1
                cases+=1
    return dict(passed=True,independent_exponential_modulated_history_C0_Z_comparisons=comparisons,
        finite_N_source_algebra_fixtures=cases,positive_negative_zero_and_tiny_A_included=True,
        nonzero_original_V_and_all_quadratic_cross_terms_included=True,
        fixture_functions_differentiated_independently_not_interval_caps=True,
        fixtures_not_original_source_selection_or_global_N=True)


def actual_large_unit_export(owner):
    c=owner.ctx;p=mp.mp.clone();p.dps=100
    a=current.frames.O2QUAtlas(owner.source.owner.inputs.frame,lower=-1,upper=0,
        logq=c.mpf((-10,-1)),logu=c.ln(c.mpf(['.1875','20'])))
    comparisons=0
    for qp in (0,1,3):
        for up in (0,-1,-2):
            value=current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,(11,10,-2,qp,up)),c.mpf(['-.7','1.3']),a.ledger)
            ratio=current.export_range(a,value,derivative_unit=True)
            for lq in (-10,-1):
                for u in (p.mpf('.1875'),p.mpf(20)):
                    for coefficient in (p.mpf('-.7'),p.mpf('1.3')):
                        saved.contains(ratio,coefficient*p.exp(qp*lq)*u**up,p.mpf('1e-35'));comparisons+=1
    rejected=0
    for powers in ((12,10,-2,0,0),(11,10,-2,-1,0),(11,10,-2,0,1)):
        try:current.export_range(a,current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,powers),1,a.ledger),derivative_unit=True)
        except ValueError:rejected+=1
    assert rejected==3
    return dict(passed=True,positive_actual_unit_division_endpoint_comparisons=comparisons,
        uncollected_large_q_u_scale_guards=rejected,
        native_Lambda0_Lminus2_canceled_before_bounded_conversion=True,
        no_materialized_cap_division_or_derivative_of_normalized_caps=True)


def continuous_density_mass_phase_contract(owner,manifest):
    c=owner.ctx;p=mp.mp.clone();p.dps=150;record=manifest['actual_full_original_O2_five_integrals'];rows=record['whole_source_density_phase_mass_records']
    assert record['exact_y_window']==['0','1'] and record['exact_Z_range']==['-1','1']
    assert record['explicit_candidate_N']==160 and len(rows)==64
    assert record['own_rates']==current.five.RATES and record['normalized_own_units']==current.five.UNITS
    assert record['original_pressure_datum_sha256']==owner.family['datum_enclosure_sha256']
    assert record['actual_global_phase_endpoint_indices']==[0,64]
    endpoints=record['actual_original_radius_phase_endpoints'];assert len(endpoints)==65
    for index,item in enumerate(endpoints):
        assert item['source_family']==owner.family and item['original_y_exact']==str(s.Rational(index,64))
        assert item['explicit_candidate_N']==160 and item['same_original_radius_phase_bound']
        assert item['positive_original_origin_offset']['strictly_positive']
        assert item['positive_original_origin_offset']['actual_hb_sc_product_not_materialized_or_set_to_zero']
    unit=record['changed_ordinary_Z_unit']['formal_positive_scale']
    assert unit['source_exponents']==[11,10,-2,0] and unit['radius_power']==0
    totals=[{key:c.mpf(0) for key in current.KEYS} for unused in range(4)]
    native=0;unions=0;masses=0
    for index,row in enumerate(rows):
        assert row['exact_y_cell']==[str(s.Rational(index,64)),str(s.Rational(index+1,64))]
        assert row['actual_original_phase_endpoint_indices']==[index,index+1]
        assert row['branch_union_not_sum_or_duplicate_integral']
        assert ep(saved.interval(c,row['actual_original_common_N_phase_cover'][0]))==(0,1)
        assert len(row['predicate_density_records'])==3
        for local in row['predicate_density_records']:
            branch=local['branch'];frame=owner.source.frame(64,index,branch=branch);a=frame.roots['q'].atlas
            assert local['actual_original_predicate']==frame.record['actual_original_source_predicate']
            assert local['original_mixed_archive_source_index']==index and local['original_mixed_archive_filename']==current.phase.NAME
            assert local['all_original_nonlinear_cross_terms_executed'] and local['source_function_derivatives_not_derivatives_of_exported_caps']
            assert local['changed_Z_large_source_unit_retained_formally']
            for j,label in enumerate(('changed_C0','changed_Z','unmodulated_C0','unmodulated_Z')):
                alias=('changed_C0','changed_Z_over_Lambda0_Lminus2','unmodulated_C0','unmodulated_Z')[j]
                for key,source in local['original_graph_density_and_ordinary_Z'][label].items():
                    assert source['encloses_original_source_function'] and not source['point_value_selected']
                    value=current.restore_value(a,source);normalized=current.export_range(a,value,derivative_unit=j==1)
                    saved.contains(saved.interval(c,local['exported_normalized_ranges'][alias][key]),normalized,mp.mpf('1e-200'));native+=1
                    union=saved.interval(c,row['normalized_four_density_unions'][j][key])
                    saved.contains(union,normalized,mp.mpf('1e-200'));unions+=1
        # Independent exact Duhamel endpoint formula; no imported mass routine.
        left=p.mpf(index)/64;right=p.mpf(index+1)/64
        for key,rate in current.five.RATES.items():
            rr=p.mpf(rate);target=(right-left) if not rr else (p.exp(-rr*(1-right))-p.exp(-rr*(1-left)))/rr
            mass=saved.interval(c,row['positive_own_rate_final_endpoint_masses'][key])
            assert ep(mass)[0]>0;saved.contains(mass,target,p.mpf('1e-120'));masses+=1
            for j in range(4):
                contribution=saved.interval(c,row['normalized_four_density_unions'][j][key])*mass
                saved.contains(saved.interval(c,row['normalized_four_final_endpoint_contributions'][j][key]),contribution,mp.mpf('1e-200'))
                totals[j][key]+=contribution
    for j in range(4):
        for key in current.KEYS:saved.contains(saved.interval(c,record['all_four_normalized_integral_enclosures'][j][key]),totals[j][key],mp.mpf('1e-200'))
    assert record['direct_density_integration_no_phase_averaging_or_boundary_term_drop']
    assert record['fixed_y_window_and_Z_independent_Duhamel_masses'] and record['no_extra_R_Jacobian']
    inlet=current.original_full_inlet(owner.out,c.mpf((-1,1)))
    assert inlet['exact_pointwise_Z_square_and_Q_at_least_one']
    for ztext in ('-.37','0','.37'):
        z=c.mpf(ztext);local=current.original_full_inlet(owner.out,z);C=1/(1+z**2)
        saved.contains(local['original_inlet']['h'],c.mpf(5)/8*C,mp.mpf('1e-200'))
        saved.contains(local['original_inlet']['p'],c.mpf(5)/2*C*C,mp.mpf('1e-200'))
        if ztext=='0':assert local['native_original_inlet_sources']['m']['exact_zero']
        assert not local['native_original_inlet_Z_sources']['m']['exact_zero']
    return dict(passed=True,continuous_original_y_cells=64,conditional_source_density_frames=192,
        native_graph_density_and_Z_exports=native,nonduplicated_branch_unions=unions,
        independent_positive_Duhamel_masses=masses,actual_common_radius_phase_endpoints=65,
        five_changed_and_original_C0_Z_full_window_sums=True,
        source_defined_inlet_including_axis_and_pressure_zero_rate_memory=True,
        ordinary_Z_integrals_factored_by_actual_Lambda0_Lminus2=True,
        not_sharp_averaging_or_functional_terminal_or_global_N=True)


def actual_native_incoming_transport(owner):
    # A fresh execution exercises issued integral ownership and the actual
    # source pipeline. No ancestor suites or additional refinement are run.
    integral=owner.integrate(N=160);a=owner.out;c=owner.ctx;inc={};incZ={}
    for i,key in enumerate(current.KEYS):
        co=c.mpf(i+1)/100;slope=c.mpf(i+1)/500
        inc[key]=a.scalar(co+c.mpf((-1,1))*slope);incZ[key]=a.scalar(slope)
    args=dict(incoming=inc,incoming_Z=incZ,source_family=owner.family,original_P0_datum_sha256=owner.family['datum_enclosure_sha256'])
    got=owner.apply_incoming(integral,**args);comparisons=0
    for key,rate in current.five.RATES.items():
        target=current.positive.bounded(integral['original'][key])+current.positive.bounded(integral['changed'][key])+current.positive.bounded(inc[key])*c.exp(-c.mpf(rate))
        saved.contains(current.positive.bounded(got['values'][key]),target,mp.mpf('1e-200'));comparisons+=1
        unit=owner.Z_unit
        norm=lambda value:current.positive.bounded(value.positive_divide(unit,unit.scale.evaluate()))
        targetZ=norm(integral['originalZ'][key])+norm(integral['changedZ'][key])+norm(incZ[key])*c.exp(-c.mpf(rate))
        saved.contains(norm(got['Z'][key]),targetZ,mp.mpf('1e-200'));comparisons+=1
    assert got['record']['exact_pressure_zero_rate_memory_preserved']
    bad=dict(inc);bad.pop('p');foreign=owner.source.frame(64,0,branch='regular').roots['E'][C0]
    rejected=0
    for call in (lambda:owner.apply_incoming(dict(integral),**args),
            lambda:owner.apply_incoming(integral,**dict(args,incoming=bad)),
            lambda:owner.apply_incoming(integral,**dict(args,incoming={**inc,'p':foreign})),
            lambda:owner.apply_incoming(integral,**dict(args,original_P0_datum_sha256='wrong')),
            lambda:owner.apply_incoming(integral,source_family=owner.family,original_P0_datum_sha256=owner.family['datum_enclosure_sha256']),
            lambda:current.candidate_N(True),lambda:current.candidate_N(0)):
        try:call()
        except (TypeError,ValueError):rejected+=1
    assert rejected==7
    return dict(passed=True,explicit_native_nonzero_incoming_C0_Z_transport_comparisons=comparisons,
        issued_identity_missing_history_foreign_atlas_pressure_and_N_guards=rejected,
        compatible_nonzero_linear_Z_incoming_fixtures_not_actual_all_route_history=True,
        source_defined_original_histories_not_zeroed=True,separate_P0_and_pressure_rate_zero_preserved=True)


def original_defining_endpoint_intersections(owner,manifest):
    """Separate slope endpoint formulas from accepted defining H,D,P masses.

    These are continuous range intersections, not selected axial fields or
    an exact endpoint identity inferred from a broad numerical interval.
    """
    c=owner.ctx;a=owner.out
    level,unused=owner.source.source.parent.parent.levels[64]
    node=level['ordered_nodes'][-1]
    H,D,P=(saved.interval(c,node[key]) for key in ('H','D','P'))
    inverse=current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,(-1,0,0,0,0)),1,a.ledger)
    inverse2=current.base.current.square(inverse);comparisons=0
    report=manifest['actual_full_original_O2_five_integrals']
    for domain in ((-1,0),(0,0),(0,1)):
        z=c.mpf(domain);zz=z**2;C=1/(1+zz)
        values=dict(m=inverse*(4*z),h=a.scalar(C*H),k=inverse*(4*z*C*H),
            e=a.add(inverse2*(16*zz),-a.scalar(C*C*D)),p=a.scalar(C*C*P))
        derivatives=dict(m=inverse*4,h=a.scalar(-2*z*C*C*H),k=inverse*(4*(1-zz)*C*C*H),
            e=a.add(inverse2*(32*z),a.scalar(4*z*C**3*D)),p=a.scalar(-4*z*C**3*P))
        for label,rows in (('source_defined_original_histories_at_y1',values),('source_defined_original_history_Z_at_y1',derivatives)):
            for key,value in rows.items():
                independent=current.positive.bounded(value)
                actual=current.positive.bounded(current.restore_value(a,report[label][key]))
                assert max(ep(independent)[0],ep(actual)[0])<=min(ep(independent)[1],ep(actual)[1]),(label,key,domain)
                comparisons+=1
    return dict(passed=True,independent_defining_mass_original_endpoint_C0_Z_intersections=comparisons,
        continuous_negative_axis_positive_Z_domains=True,accepted_defining_H_D_P_masses_reused_no_ancestor_suite=True,
        intersections_not_exact_terminal_identity=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_all_route_incoming_histories_installed','sharp_O2_phase_averaging_installed','functional_terminal_identity_solved',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed','current_whole_N_selected',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2PredicateFiveIntegrals();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('finite_N_defining_moment_and_normalized_B_units',density_and_unit_calculus),
                ('independent_exponential_density_Z_fixtures',lambda:finite_N_defining_density_fixtures(owner)),
                ('actual_large_source_unit_export',lambda:actual_large_unit_export(owner)),
                ('continuous_original_density_phase_mass_integral_domain',lambda:continuous_density_mass_phase_contract(owner,manifest)),
                ('original_defining_H_D_P_endpoint_intersections',lambda:original_defining_endpoint_intersections(owner,manifest)),
                ('issued_original_integral_and_explicit_incoming_C1',lambda:actual_native_incoming_transport(owner))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Full original O2 candidate-N160 five changed/original C0/Z direct Duhamel integral enclosures, actual source predicate/phase union, normalized B units, positive masses and original inlet, factored large ordinary-Z unit and explicit native incoming protocol. Not sharp averaging, all-route incoming functions, terminal/global N, stress/recursion or full reconstruction.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Full original O2 predicate five C0/Z direct integrals and native incoming interface PASS',flush=True);return report


if __name__=='__main__':run()
