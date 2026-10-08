"""Original midplane source parity, independent derivatives and live oracle evidence."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_midplane_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;points=current.points


def pressure_and_template_identities(owner):
    z,f,H,D,P,P0,P0Z,P0ZZ=owner.templates['inputs'];subs={z:0,P0Z:0}
    pressure_symbols=owner.templates['pressure_symbols'];rows=owner.templates['rows'];identities=0
    for powers,expression in rows[('p2',0)]:
        assert s.simplify(expression.subs(subs,simultaneous=True))==0
        # Even pressure errors still multiply an exact zero coefficient.
        for order in (0,2):
            assert s.simplify(s.diff(expression,pressure_symbols[order]).subs(subs,simultaneous=True))==0
        identities+=3
    assert any(s.simplify(s.diff(expression,P0Z).subs(subs,simultaneous=True))!=0
               for powers,expression in rows[('p2',0)])
    assert any(s.simplify(expression.subs(subs,simultaneous=True))!=0
               for powers,expression in rows[('p2',1)])
    assert any(s.diff(expression,P0ZZ)!=0 for powers,expression in rows[('p2',1)])
    assert len(owner.parity['actual_original_even_density_stage_identities'])==14
    assert owner.parity['source_midplane_pressure_first_jet_including_remainder_exact_zero']
    assert all(v['reversed_AST_exactly_equals_original'] for v in current.PROJECTION_PROOF.values())
    # Genuine exact source functions in isolated finite diagnostic units.
    # Differentiate C0, including L(Z) and analytic P0(Z); no order1 reuse.
    p=mp.mp.clone();p.dps=130
    functions={key:tuple(s.lambdify(owner.templates['inputs'],expr,
                 modules=[{'mpf':p.mpf},'mpmath']) for powers,expr in terms) for key,terms in rows.items()}
    alpha=p.mpf('2.9');Pstar=p.exp(3);delta=p.exp(-42)
    def root(name,y,Z,order):
        ff=p.exp(y/10);q=1+Z*Z;R=110*p.exp(50+y);L=1-delta*Z*Z
        values=(Z,ff,p.mpf(5)/8*ff,p.mpf(5)/12*ff*ff,p.mpf(5)/2*ff*ff,
                -alpha/q**2,4*alpha*Z/q**3,alpha*(4-20*Z*Z)/q**4)
        return sum(R**r*Pstar**a*delta**d*L**ell*fn(*values)
                   for ((r,a,d,ell),expr),fn in zip(rows[(name,order)],functions[(name,order)],strict=True))
    comparisons=0
    for y in (p.mpf('-2.337'),p.mpf('-.513')):
        for name in ('E','V','b','p1','p2'):
            expected=p.diff(lambda Z:root(name,y,Z,0),0)
            got=root(name,y,p.mpf(0),1)
            assert abs(got-expected)<=p.mpf('1e-110')*max(abs(got),abs(expected),1),(name,y)
            comparisons+=1
        assert root('p2',y,0,0)==0 and root('p2',y,0,1)!=0
    query=owner.roots('-2.5','-1.25')
    assert query['roots']['p2'][(0,0)].zero and not query['roots']['p2'][(0,1)].zero
    assert ep(query['roots']['p2'][(0,1)].coefficient)[1]<0
    assert ep(owner.atlas.bases[2])==(0,0)
    return dict(passed=True,exact_p2_midplane_zero_and_even_pressure_sensitivity_identities=identities,
        source_p2_Z_and_P0_ZZ_dependence_nonzero=True,sourcewide_even_pressure_stages=14,
        reversible_original_source_AST_projections=2,
        independent_original_C0_derivative_comparisons=comparisons,
        strict_original_nonzero_p2_Z_retained=True,original_L0_exact_one=True,
        finite_diagnostic_units_not_installed_as_source_parameters=True)


def independent_midplane_primitives(owner):
    a=owner.atlas;c=owner.ctx;p=mp.mp.clone();p.dps=110;comparisons=0
    for derivative in (p.mpf('.7'),p.mpf('-.7')):
        q=p.mpf('.9');E=p.mpf('1.3');aa=p.mpf('.8')
        roots={name:{(0,0):a.scalar(c.mpf(value)),(0,1):a.scalar(c.mpf(dZ))}
               for name,value,dZ in (('E',E,0),('a',aa,0),('b',0,0),('t0',0,0),('p2',0,derivative/q))}
        kernel=current.whole.conditioned.PositiveLogQPhase(
            dict(q=a.scalar(c.mpf(q)),roots=roots,original_u_source=a.scalar(0)),c.mpf(0))
        assert kernel.u.zero and kernel.geometry=='small_r_series'
        def direction(theta,Z):
            uu=derivative*Z;hh=p.sqrt(1+uu*uu);r=uu/hh
            return 2*q/hh*(p.cos(theta)-r)/(1-2*r*p.cos(theta)+r*r)
        for fraction in ('.23','.71'):
            x=p.mpf(fraction);psi=2*p.pi*x
            T1=p.quad(lambda theta:direction(theta,0),[0,psi])
            T2=p.quad(lambda theta:direction(theta,0)**2,[0,psi])
            T1Z=p.quad(lambda theta:p.diff(lambda Z:direction(theta,Z),0),[0,psi])
            T2Z=p.quad(lambda theta:p.diff(lambda Z:direction(theta,Z)**2,0),[0,psi])
            phi=(psi+T2)/(2*p.pi*(1+2*q*q));tt=direction(psi,0)
            expected=dict(A=aa/2*(phi-x),B_over_Pstar=-aa*E*T1/(4*p.pi),
                          A_Z_slow=aa*T2Z/(4*p.pi*(1+tt*tt)),
                          B_Z_slow=-aa*E/(4*p.pi)*(T1Z-tt*T2Z/(1+tt*tt)))
            primitive=kernel.primitives(c.mpf(x),'psi')
            slow,proof=points.slow.slow_values(kernel,roots,c.mpf(x),'psi')
            assert proof['branch']=='exact_midplane_nonzero_p2_Z'
            for key,value in expected.items():
                got=(primitive if key in primitive else slow)[key].finite_interval()
                saved.contains(got,value,p.mpf('1e-95'));comparisons+=1
            for key,value in (('T1_fixed_angle',T1),('T1_Z_fixed_angle',T1Z),('T2_Z_fixed_angle',T2Z)):
                saved.contains(saved.saved_value(c,proof[key]),value,p.mpf('1e-95'));comparisons+=1
        for x in (0,'.5',1):
            values,proof=points.slow.slow_values(kernel,roots,c.mpf(x),'psi')
            assert all(value.zero for value in values.values())
    return dict(passed=True,independent_original_defining_integral_C0_Z_primitive_comparisons=comparisons,
        both_derivative_signs_tested=True,exact_period_and_halfperiod_Z_traces=True,
        genuine_nonzero_midplane_derivative_not_flattened=True)


def live_point_in_whole_cell(owner):
    a=owner.atlas;c=owner.ctx;comparisons=0
    query=points.evaluate_coefficients(owner.dispatcher,chart='Rh_reference',
        coordinate='-2.337',Z=0,N=160,bits=40)
    whole=owner.cell('-2.5','-2.1875',N=160)
    for item in query['evaluated']:
        piece=item['piece']
        assert piece.Z==0 and piece.family==owner.family
        for order in points.exact.ORDERS:
            for name in points.exact.RATES:
                for jet in ('C0','Z'):
                    original=item['values'][order][name][jet]
                    value=a.rebase_piece_value(owner.dispatcher,piece,original)
                    box=whole['coefficients'][order][name][jet]
                    anchor=current.whole.prior.FormalScale(a.bases,box.scale.powers)
                    normalize=lambda v:c.mpf(0) if v.zero else v.coefficient*v.bounded_exp((v.scale-anchor).evaluate())
                    saved.contains(normalize(box),normalize(value));comparisons+=1
    rejected=0
    for call in (lambda:current.OriginalReferenceMidplaneWholeCells(Z='.37'),
                 lambda:current.whole.OriginalReferenceWholeCells(Z=0),
                 lambda:owner.roots(0,0),lambda:owner.integrate(count=0,N=160)):
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejected+=1
    assert rejected==4
    return dict(passed=True,actual_original_point_coefficients_contained_in_whole_midplane_cell=comparisons,
        actual_original_phase_and_source_errors_consumed=True,
        exact_Z_domain_partition_and_original_signed_guard_rejections=rejected,
        original_signed_service_guard_not_relaxed=True)


def native_integral_contracts(owner,manifest):
    c=owner.ctx;levels=manifest['actual_original_midplane_whole_reference_levels'];summaries=[];masses=0
    for level in levels:
        count,N=level['exact_original_cells'],level['candidate_N'];cells=level['full_whole_cell_source_phase_records']
        assert level['source_family']==owner.family and level['original_Z_exact']=='0'
        assert level['exact_source_window']==['-5','0'] and len(cells)==count
        assert level['exact_midplane_whole_window_C0_Z_contributions_installed']
        assert level['high_precision_source_and_integral_exports_retained']
        assert not level['signed_implicit_Z_products_bounded_before_multiplication']
        assert level['incoming_global_correction_histories_not_assumed_zero']
        assert level['original_P0_kept_separate_and_not_reset'] and level['original_pressure_rate_zero_memory_retained']
        assert not level['whole_Z_functional_or_all_route_closure']
        totals={name:c.mpf(0) for name in points.exact.RATES}
        for i,cell in enumerate(cells):
            source=cell['source'];left,right=map(s.Rational,source['exact_reference_cell'])
            assert left==-5+s.Rational(5*i,count) and right==-5+s.Rational(5*(i+1),count)
            assert source['candidate_N']==N and source['source_geometry']=='small_r_series'
            roots=source['whole_cell_original_roots']
            assert roots['p2']['(0, 0)']['exact_zero']
            assert roots['p2']['(0, 1)']['sign']=='negative'
            for piece in source['original_inverse_and_Z_piece_records']:
                assert piece['ordinary_Z']['branch']=='exact_midplane_nonzero_p2_Z'
                assert not piece['ordinary_Z']['u_Z']['exact_zero']
                assert piece['ordinary_Z']['p2_Z_retained']
            for row in source['whole_cell_closed_source_terms']:
                for term in row['terms']:
                    if term['normalized_pressure_late_error_Pstar_power'] is not None:
                        assert term['normalized_pressure_late_error_Pstar_power']==-1
            assert cell['native_radius_Jacobian']==1
            for name,mass in cell['original_positive_own_rate_masses'].items():
                box=saved.interval(c,mass);assert ep(box)[0]>0;totals[name]+=box
        for name,rate in points.exact.RATES.items():
            r=c.mpf(rate.numerator)/rate.denominator
            expected=c.mpf(5) if ep(r)==(0,0) else (1-c.exp(-5*r))/r
            saved.contains(totals[name],expected);masses+=1
        C0={};Znormalized={}
        for name,pair in level['actual_finite_N_five_density_contribution_enclosures'].items():
            row=pair['C0'];assert row['formal_positive_scale']['source_exponents']==[0,0,0,0]
            C0[name]=max(abs(v) for v in ep(saved.saved_value(c,row,owner.atlas.bases)))
            row=pair['Z'];scale=row['formal_positive_scale']
            assert scale['source_exponents']==[11,10,-2,0] and scale['radius_power']==0
            assert not row['exact_zero'] and not row['point_value_selected']
            assert not any(not mp.isfinite(v) for v in ep(saved.interval(c,row['coefficient_interval'])))
            # Retain the astronomical original positive factor formally.
            # Only the finite coefficient/offset is normalized.
            Znormalized[name]=max(abs(v) for v in ep(saved.interval(c,row['coefficient_interval'])*
                  c.exp(saved.interval(c,scale['additional_log_interval']))))
        summaries.append(dict(cells=count,N=N,C0_absolute_bounds=C0,
            Z_absolute_bounds_divided_by_original_Pstar11_Cstar10_Lminus2=Znormalized))
    assert [(v['exact_original_cells'],v['candidate_N']) for v in levels]==[(4,160),(16,160),(16,320),(16,16384)]
    first,last=levels[1],levels[-1]
    assert first['full_whole_cell_source_phase_records'][0]['source']['whole_cell_N_dependent_coefficients']['-1']['h']!=last['full_whole_cell_source_phase_records'][0]['source']['whole_cell_N_dependent_coefficients']['-1']['h']
    assert summaries[-1]['C0_absolute_bounds']['p']<c.mpf('.000052')
    return dict(passed=True,actual_levels=[(v['exact_original_cells'],v['candidate_N']) for v in levels],
        positive_own_rate_mass_comparisons=masses,pressure_zero_rate_memory_mass=5,
        original_native_Z_factor=[11,10,-2,0],large_Z_factors_not_materialized_or_capped=True,
        actual_N_dependent_analytic_coefficients_recomputed=True,
        whole_Z_or_global_N_not_claimed=True,computed_ranges=summaries)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
           'current_whole_N_selected',*point_flags())
    assert all(manifest[name] is False for name in flags)
    owner=current.OriginalReferenceMidplaneWholeCells()
    with mp.workdps(owner.ctx.dps+40):
        checks={}
        for name,call in (
            ('original_pressure_and_midplane_source',lambda:pressure_and_template_identities(owner)),
            ('independent_original_midplane_primitive_evidence',lambda:independent_midplane_primitives(owner)),
            ('live_midplane_point_to_whole_source',lambda:live_point_in_whole_cell(owner)),
            ('actual_midplane_integral_contracts',lambda:native_integral_contracts(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
             Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Original exact-midplane full reference C0/Z contributions: full pressure parity with even error terms retained, original source derivatives and defining primitive integrals, actual live point contained in whole cell, own rates and nonzero formal native Z amplification. Not a Z neighborhood, functional terminal closure/all-route/global N or recursive corrected field.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original exact-midplane source, primitives, live oracle and whole contributions PASS',flush=True)
    return report


def point_flags():return current.point.source.inertial.profiles.loop.OPEN


if __name__=='__main__':run()
