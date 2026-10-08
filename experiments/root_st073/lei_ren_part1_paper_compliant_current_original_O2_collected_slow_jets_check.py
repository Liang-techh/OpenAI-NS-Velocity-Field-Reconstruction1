"""Direct carrier, signed monomial bounds and actual native slow norm check."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_collected_slow_jets as current
import lei_ren_part1_paper_compliant_current_original_O2_uniform_mean_bias_check as previous

saved=previous.saved;ep=current.ep;C0,Y,Z,YZ=current.mixed.ORDERS


def exact_calculus_and_finite_signed_densities(owner):
    c=owner.ctx;frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas
    p=mp.mp.clone();p.dps=150;y0,z0=p.mpf('.41'),p.mpf('-.37');comparisons=0;cases=0
    E,V,A,B,t=s.symbols('E V A B t',real=True)
    def density(E,V):return dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
    changed=density(E*s.exp(t*A),V+t*B);old=density(E,V)
    target=dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)
    for key in current.integrals.KEYS:assert s.simplify(s.diff(changed[key]-old[key],t).subs(t,0)-target[key])==0
    for vv in ('-.21','.21'):
        for aa in ('-.5','0','1e-50','.5'):
            data=dict(E=('1.1','.2','-.3','.04'),V=(vv,'.07','.4','-.05'),
                A=(aa,'.17','-.13','.11'),B=('.53','-.09','.23','-.08'))
            jets={name:current.MixedJet(a,{order:a.scalar(c.mpf(value)) for order,value in zip(current.mixed.ORDERS,rows,strict=True)})
                for name,rows in data.items()}
            actual=current.leading_densities(a,**jets)
            def f(name,y,z):
                v,dy,dz,dyz=map(p.mpf,data[name]);return v+dy*(y-y0)+dz*(z-z0)+dyz*(y-y0)*(z-z0)
            def first(key,y,z):
                E,V,A,B=(f(name,y,z) for name in ('E','V','A','B'))
                return dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)[key]
            for key,jet in actual.items():
                values={C0:first(key,y0,z0),Y:p.diff(lambda y:first(key,y,z0),y0),
                    Z:p.diff(lambda z:first(key,y0,z),z0),YZ:p.diff(lambda y:p.diff(lambda z:first(key,y,z),z0),y0)}
                for order,value in values.items():saved.contains(jet[order].finite_interval(),value,p.mpf('1e-120'));comparisons+=1
            cases+=1
    assert current.source_unit_theorem()['exact_direct_original_u_y_Z_yZ_identities']==3
    return dict(passed=True,exact_signed_first_order_exponential_linearization_identities=5,
        independent_compatible_first_order_function_cases=cases,independent_four_slow_density_derivative_comparisons=comparisons,
        exact_direct_source_u_product_identities=3,nonzero_V_and_signed_energy_terms_retained=True)


def finite_signed_monomial_and_unit_bounds(owner):
    """Independent finite diagnostic atlas tests reciprocal inequalities."""
    c=owner.ctx;p=mp.mp.clone();p.dps=150;checks=0;unitchecks=0;guardchecks=0
    actual=owner.source.frame(64,0,branch='regular');ledger=dict(actual.roots['q'].atlas.ledger)
    P,C,L,d,G=map(p.mpf,('1.01','1.02','.94','.1','2'));Lambda=P**11*C**10;guard=p.mpf(3)/16
    powers=((0,0,0,1,-5),(0,0,0,-2,-1),(-1,0,1,-.5,-.5),
        (0,0,0,-1,1),(0,0,0,1,2),(1,-1,-1,0,0))
    for branch in ('positive','negative'):
        for q in map(p.mpf,('.03','.6')):
            for u in (guard,Lambda*G*q/d*p.mpf('.6')):
                assert guard<=u<=Lambda*G*q/d
                bases=tuple(c.ln(c.mpf(str(v))) for v in (P,C,L,q,u))
                a=SimpleNamespace(ctx=c,bases=bases,ledger=ledger)
                a.scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
                dd=current.prior.ScaledEnclosure(current.prior.FormalScale(bases,offset=c.ln(c.mpf(str(d)))),1,ledger)
                frame=SimpleNamespace(branch=branch,roots={'q':SimpleNamespace(atlas=a)},kernel=SimpleNamespace(dstar=dd))
                def describe(value):
                    if value is not frame:raise ValueError('Finite fixture owner required')
                    return {}
                fixture=SimpleNamespace(source=SimpleNamespace(describe=describe,g=c.mpf((-2,-2))))
                for row in powers:
                    value=current.prior.ScaledEnclosure(current.prior.FormalScale(bases,row,offset=c.mpf('-.2')),c.mpf('-.73'),ledger)
                    bound=current.absolute_monomial_bound(fixture,frame,value)
                    exact=p.mpf('.73')*p.exp(p.mpf('-.2'))
                    for v,power in zip((P,C,L,q,u),row,strict=True):exact*=v**power
                    upper=ep(bound.finite_interval())[1]
                    saved.contains(c.mpf((0,upper)),exact,p.mpf('1e-110'));checks+=1
                for order,unit in current.UNITS.items():
                    row=unit[:3]+(1,-1)
                    value=current.prior.ScaledEnclosure(current.prior.FormalScale(bases,row),c.mpf('-.73'),ledger)
                    norm=current.native_norm(fixture,frame,value,order)
                    saved.contains(norm['cap'],p.mpf('.73')*q/u,p.mpf('1e-110'));unitchecks+=1
                frame.branch='regular'
                for row in ((0,0,0,-1,0),(0,0,0,0,-1)):
                    value=current.prior.ScaledEnclosure(current.prior.FormalScale(bases,row),1,ledger)
                    try:current.absolute_monomial_bound(fixture,frame,value)
                    except ValueError:guardchecks+=1
    assert checks==48 and unitchecks==32 and guardchecks==16
    return dict(passed=True,independent_signed_absolute_monomial_inequality_checks=checks,
        finite_common_unit_exports=unitchecks,regular_reciprocal_q_or_u_rejections=guardchecks,
        boundary_abs_u_and_negative_half_integer_powers_tested=True,
        finite_fixture_atlas_not_original_source_replacement=True)


def actual_source_binding_and_norms(owner,manifest):
    c=owner.ctx;rows=manifest['whole_original_O2_collected_slow_source_cells']
    old=json.loads(gzip.decompress((current.HERE/current.mixed.NAME).read_bytes()))['whole_original_O2_correlated_qy_mixed_density_cells']
    assert manifest['source_family']==owner.family and manifest['original_source_cells']==64
    assert manifest['original_predicate_frames']==192 and len(rows)==64
    assert manifest['exact_outer_source_domain']==dict(y=['0','1'],Z=['-1','1'],phi=['0','1'])
    counts=dict(frames=0,raw=0,unchanged_phase=0,source=0,norms=0,direct=0)
    maxima={name:{str(order):c.mpf(0) for order in current.mixed.ORDERS} for name in ('A','B',*current.integrals.KEYS)}
    for index,row in enumerate(rows):
        assert set(row)==set(current.integrals.BRANCHES)
        for branch,item in row.items():
            frame=owner.source.frame(64,index,branch=branch);a=frame.roots['q'].atlas
            assert frame.roots['q'][C0] is frame.kernel.q
            assert item['source']['original_pressure_q_eta_nu_and_P0_unchanged']
            assert item['source']['independently_selected_compatible_C0_derivative_fields_claimed'] is False
            assert item['source']['regular_raw_p2_replaces_dstar_u_over_q_for_derivatives']==(branch=='regular')
            assert item['source']['signed_source_kernel_and_u_unchanged']==(branch!='regular')
            raw=owner.source.source.source_frame(64,index,Z_lower=a.bounds[0],Z_upper=a.bounds[1])
            for order,value in raw.roots['p2'].rows.items():
                rebased=current.frames.rebase_original_root(a,value,raw_frame=raw,source_owner=owner.source.source)
                previous.exact_native(c,item['source']['raw_original_p2_C0_y_Z_yZ'][str(order)],current.base.encoded(rebased.record()));counts['raw']+=1
                if branch=='regular':previous.exact_native(c,frame.roots['p2'][order].record(),rebased.record())
            if branch=='regular':
                direct=current.direct_u_jets(a,frame.roots['p2'],frame.roots['q'],frame.kernel.dstar,frame.u[C0])
                for order in (Y,Z,YZ):previous.exact_native(c,frame.u[order].record(),direct[order].record());counts['direct']+=1
                assert max(abs(v) for v in ep(frame.u[C0].finite_interval()))<=c.mpf(1)/4
                assert frame.u[Y].scale.powers[3]>=0 and frame.u[YZ].scale.powers[3]>=0
            replay=owner.collect(frame);again=current.base.encoded(replay['record'])
            unstrengthened=item['actual_full_phase_mixed_before_C0_strengthening']['complete_original_A_B_mixed']
            predecessor=old[index][branch]['actual_whole_phase_mixed']['complete_original_A_B_mixed']
            for name in ('A','B'):
                for order in (C0,Z):
                    previous.exact_native(c,unstrengthened[name][str(order)],predecessor[name][str(order)]);counts['unchanged_phase']+=1
            for name in ('actual_same_source_A_B_C0_y_Z_yZ','native_same_source_signed_first_order_density_C0_y_Z_yZ'):
                for component,orders in item[name].items():
                    for order,value in orders.items():previous.exact_native(c,value,again[name][component][order]);counts['source']+=1
            for name in ('primitive_native_norms','first_order_density_native_norms'):
                for component,orders in item[name].items():
                    for order,record in orders.items():
                        other=again[name][component][order]
                        assert record['bound_only_not_source_function_or_derivative_of_cap']
                        assert record['positive_envelope_lower_endpoint_not_source_magnitude_lower_bound']
                        previous.exact_interval(c,record['normalized_absolute_cap'],other['normalized_absolute_cap'])
                        for key in ('collected_absolute_native_bound','exact_positive_normalization_unit'):previous.exact_native(c,record[key],other[key])
                        unit=record['exact_positive_normalization_unit']['formal_positive_scale']
                        assert tuple(unit['source_exponents'])+(unit['radius_power'],)=={str(k):v for k,v in current.UNITS.items()}[order]
                        upper=max(ep(maxima[component][order])[1],ep(saved.interval(c,record['normalized_absolute_cap']))[1])
                        maxima[component][order]=c.mpf((0,upper));counts['norms']+=1
            assert item['exact_leading_density_true_phase_mean_and_slow_derivatives_zero']
            assert item['same_source_norms_not_point_values'] and item['original_radius_common_N_phase_not_yet_spatially_integrated']
            counts['frames']+=1
        if (index+1)%16==0:print('Checked original native slow norms:',index+1,'/64',flush=True)
    for name,orders in maxima.items():
        for order,value in orders.items():previous.exact_interval(c,value,manifest['whole_source_normalized_slow_norms'][name][order])
    copied=replace(owner.source.frame(64,0,branch='regular'))
    try:owner.source.describe(copied)
    except ValueError:counts['copied_owner_rejected']=True
    assert counts.get('copied_owner_rejected')
    return dict(passed=True,**counts,original_P0_and_parameter_source_unchanged=True,
        all_native_raw_p2_source_rows_bound_to_accepted_source=True,
        dominating_normalization_units_not_exact_raw_derivative_powers=True,
        continuous_named_source_predicates_not_selected_compatible_fields=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_zero_mean_phase_primitive_evaluator_installed','oscillatory_spatial_integral_remainder_enclosed',
        'sharp_O2_phase_averaging_installed','actual_all_route_incoming_histories_installed','functional_terminal_identity_solved',
        'current_whole_N_selected','all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2CollectedSlowJets();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('independent_direct_and_signed_first_order_calculus',lambda:exact_calculus_and_finite_signed_densities(owner)),
                ('independent_signed_native_monomial_and_unit_bounds',lambda:finite_signed_monomial_and_unit_bounds(owner)),
                ('actual_source_carrier_phase_and_slow_norm_binding',lambda:actual_source_binding_and_norms(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original O2 raw regular p2 direct source derivatives, unchanged phase C0/Z, all actual first-order five density slow jets and finite native dominating norms. Not selected interval field values, evaluated phase primitives, actual spatial remainder, terminal/global N or scale recursion.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Original O2 collected native slow jets PASS',flush=True);return report


if __name__=='__main__':run()
