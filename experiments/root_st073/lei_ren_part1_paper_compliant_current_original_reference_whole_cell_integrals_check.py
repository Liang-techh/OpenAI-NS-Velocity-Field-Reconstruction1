"""Focused factor, original whole-cell, implicit-product and own-rate evidence."""
from dataclasses import replace
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals as current

points=current.points;base=current.base;ep=current.ep


def interval(c,row):
    if hasattr(row,'_mpi_'):return c.mpf(ep(row))
    return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                  mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))


def saved_value(c,row,bases=None):
    scale=row['formal_positive_scale']
    powers=(*scale['source_exponents'],scale['radius_power'])
    if any(powers):assert bases is not None
    log=interval(c,scale['additional_log_interval'])
    if bases is not None:log+=sum((b*p for b,p in zip(bases,powers) if p),c.mpf(0))
    return interval(c,row['coefficient_interval'])*c.exp(log)


def contains(box,value,tolerance=0):
    lo,hi=ep(box)
    if hasattr(value,'_mpi_'):
        a,b=ep(value);assert lo-tolerance<=a<=b<=hi+tolerance,(lo,hi,a,b)
    else:assert lo-tolerance<=value<=hi+tolerance,(lo,hi,value)


def exact_factor_and_ordering(owner):
    atlas=owner.atlas;c=atlas.ctx
    P,C,L,y,o=s.symbols('logP logC logL y offset',real=True)
    cases=((0,0,0,0,0),(1,-2,1,0,1),(-3,2,-1,0,-2),
           (2,-1,.5,0,1),(-1,0,-.5,0,0))
    p=mp.mp.clone();p.dps=150;diagnostics=0
    for powers in cases:
        a,d,ell,_,r=powers
        original=a*P+d*(-4*P-30)+ell*L+r*(s.log(110)+10*C+10*P+y)+o
        collected=(a-4*d+10*r)*P+10*r*C+ell*L+o-30*d+r*(s.log(110)+y)
        assert s.expand(original-collected)==0
        mapped=atlas.scale(powers,c.mpf(['-.7','-.2']),c.mpf('.3'))
        assert mapped.powers==(a-4*d+10*r,10*r,ell,0,0)
        # Finite units test the affine map independently. They never enter
        # an original source frame, manifest, frequency or field evaluation.
        logP,logC,z,yy,offset=p.mpf(3),p.mpf(2),p.mpf('.37'),p.mpf('-.4'),p.mpf('.3')
        logdelta=-4*logP-30;logL=p.log(1-p.exp(logdelta)*z*z)
        logR=p.log(110)+10*(logC+logP)+yy
        before=p.exp(a*logP+d*logdelta+ell*logL+r*logR+offset)
        after=p.exp(mapped.powers[0]*logP+mapped.powers[1]*logC+ell*logL+
                    offset-30*d+r*(p.log(110)+yy))
        assert abs(before-after)<=p.mpf('1e-140')*max(abs(before),abs(after),1)
        diagnostics+=1
    # The actual selected logC common factor is astronomically larger than
    # logP. The source difference must be collected before scale ordering.
    large=atlas.term((1,0,0,0),1,coordinate=0)
    small=atlas.term((1,0,1,0),1,coordinate=0)
    got=atlas.add(small,-large)
    normalized=got.coefficient*got.bounded_exp((got.scale-large.scale).evaluate())
    delta_cover=small.bounded_exp((small.scale-large.scale).evaluate())
    assert ep(normalized)[1]<0
    contains(normalized,-1+delta_cover)
    assert ep(delta_cover)[0]==0 and ep(delta_cover)[1]>0
    assert atlas.ledger['collected_relative_scale_additions']>0
    assert atlas.ledger['directed_independent_log_rescalings']==0
    assert ep(atlas.bases[2])[0]<=0<=ep(atlas.bases[2])[1]
    return dict(passed=True,exact_affine_power_identities=len(cases),
        independent_finite_unit_diagnostics=diagnostics,
        actual_huge_shared_radius_delta_ordering_and_positive_tail_regression=True,
        no_native_radius_exponential_or_huge_log_subtraction=True)


def implicit_product_references(owner):
    t=s.Symbol('t',real=True)
    assert s.cancel(s.Rational(1,2)-t/(1+t*t)-(t-1)**2/(2*(1+t*t)))==0
    assert s.cancel(s.Rational(1,2)+t/(1+t*t)-(t+1)**2/(2*(1+t*t)))==0
    a=owner.atlas;c=a.ctx;p=mp.mp.clone();p.dps=110
    comparisons=0
    for sign in (1,-1):
        roots={name:{(0,0):a.scalar(value),(0,1):a.scalar(derivative)}
               for name,value,derivative in (('E',c.mpf('1.3'),c.mpf('-.2')),
                   ('a',c.mpf('.8'),0),('b',0,0),('t0',0,0),('p2',2*sign,c.mpf('.7')))}
        q=a.scalar(c.mpf('.9'));u=roots['p2'][(0,0)]*q
        kernel=current.conditioned.PositiveLogQPhase(dict(q=q,roots=roots,original_u_source=u),c.mpf(0))
        assert kernel.geometry=='signed_Mobius'
        # Independent defining integrals T1=int(t), T2=int(t^2) in finite
        # diagnostic units; no transformed T1/T2 source expression is reused.
        def direction(theta,Z):
            uu=(2*sign+p.mpf('.7')*Z)*p.mpf('.9');hh=p.sqrt(1+uu*uu);rr=uu/hh
            return 2*p.mpf('.9')/hh*(p.cos(theta)-rr)/(1-2*rr*p.cos(theta)+rr*rr)
        for x in ('.23','.71'):
            uu=2*sign*p.mpf('.9');rr=uu/p.sqrt(1+uu*uu);chi=2*p.pi*p.mpf(x)
            psi=2*p.atan2((1-rr)*p.sin(chi/2),(1+rr)*p.cos(chi/2))
            if psi<0:psi+=2*p.pi
            T1=p.quad(lambda theta:direction(theta,0),[0,psi])
            T1Z=p.diff(lambda Z:p.quad(lambda theta:direction(theta,Z),[0,psi]),0)
            T2Z=p.diff(lambda Z:p.quad(lambda theta:direction(theta,Z)**2,[0,psi]),0)
            bounded,proof=current.bounded_signed_implicit_B_Z(a,kernel,roots,c.mpf(x),'E')
            for key,expected in (('T1_fixed_angle',T1),('T1_Z_fixed_angle',T1Z),('T2_Z_fixed_angle',T2Z)):
                contains(saved_value(c,proof[key]),expected,p.mpf('1e-95'));comparisons+=1
            tt=direction(psi,0)
            expected=-p.mpf('.8')/(4*p.pi)*(-p.mpf('.2')*T1+
                p.mpf('1.3')*(T1Z-tt*T2Z/(1+tt*tt)))
            contains(bounded.finite_interval(),expected,p.mpf('1e-95'));comparisons+=1
            assert not proof['selected_endpoint_midpoint_or_source_cap_value']
    return dict(passed=True,exact_signed_rational_factor_certificates=2,
        independent_defining_integral_C0_Z_and_fixed_phase_B_Z_comparisons=comparisons,
        positive_and_negative_signed_Mobius_cases=True,
        finite_units_are_diagnostics_only_not_original_source_values=True)


def live_rebase_and_whole_cell(owner):
    a=owner.atlas;dispatcher=owner.dispatcher;rebases=compared=0;last=None
    with mp.workdps(a.ctx.dps+40):
        for chart,coordinate in (('Rh_reference','-2.337'),('O2_slope','.53')):
            query=points.evaluate_coefficients(dispatcher,chart=chart,coordinate=coordinate,Z='.37',N=160)
            for item in query['evaluated']:
                piece=item['piece']
                for value in piece.values.values():
                    got=a.rebase_piece_value(dispatcher,piece,value)
                    assert got.scale.bases is a.bases and got.ledger is a.ledger
                    rebases+=1
                if chart=='Rh_reference':
                    whole=owner.cell('-2.5','-1.25',N=160)
                    for order in points.exact.ORDERS:
                        for name in points.exact.RATES:
                            for jet in ('C0','Z'):
                                box=whole['coefficients'][order][name][jet]
                                value=a.rebase_piece_value(dispatcher,piece,item['values'][order][name][jet])
                                anchor=current.prior.FormalScale(a.bases,box.scale.powers)
                                def normalized(v):
                                    if v.zero:return a.ctx.mpf(0)
                                    return v.coefficient*v.bounded_exp((v.scale-anchor).evaluate())
                                contains(normalized(box),normalized(value));compared+=1
                last=piece
    rejected=0
    other=current.factors.OriginalSourceFactorAtlas(a.frame,Z='.38')
    bad_ledger=current.prior.ScaledEnclosure(last.values['all_N_original_E_C0'].scale,1,{})
    calls=(lambda:a.rebase_piece_value(dispatcher,replace(last),last.values['all_N_original_E_C0']),
           lambda:other.rebase_piece_value(dispatcher,last,last.values['all_N_original_E_C0']),
           lambda:a.rebase_piece_value(dispatcher,last,bad_ledger),
           lambda:a.parameter('Cstar'),lambda:current.OriginalReferenceWholeCells(Z=0),
           lambda:owner.roots(0,0),lambda:owner.integrate(count=0,N=160))
    for call in calls:
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejected+=1
    assert rejected==len(calls)
    return dict(passed=True,issued_original_two_chart_factor_rebases=rebases,
        actual_point_C0_Z_coefficients_contained_in_whole_reference_cell=compared,
        unissued_frame_Z_ledger_parameter_midplane_and_partition_rejections=rejected,
        closed_whole_cell_formulas_not_point_extrapolation=True)


def whole_window_contracts(owner,manifest):
    c=owner.ctx;p=mp.mp.clone();p.dps=110;mass_checks=0;widths=[]
    Z=s.Symbol('Z',real=True);alpha=s.Symbol('alpha',positive=True)
    P0=-alpha/(1+Z*Z)**2
    assert s.simplify(s.diff(P0,Z)-4*alpha*Z/(1+Z*Z)**3)==0
    assert s.simplify(s.diff(P0,Z,2)-alpha*(4-20*Z*Z)/(1+Z*Z)**4)==0
    linear=0
    for terms in owner.templates['rows'].values():
        for _,expr in terms:
            for P in owner.templates['pressure_symbols']:
                for Q in owner.templates['pressure_symbols']:assert s.diff(expr,P,Q)==0
            linear+=1
    for level in manifest['original_reference_whole_window_refinements']:
        count=level['exact_original_cells'];cells=level['full_whole_cell_source_phase_records']
        assert len(cells)==count and level['exact_source_window']==['-5','0']
        totals={k:c.mpf(0) for k in points.exact.RATES};late_terms=0
        for i,cell in enumerate(cells):
            source=cell['source']
            left,right=map(s.Rational,source['exact_reference_cell'])
            assert left==-5+s.Rational(5*i,count) and right==-5+s.Rational(5*(i+1),count)
            assert source['source_family']==owner.family and source['candidate_N']==160
            assert cell['native_radius_Jacobian']==1 and source['source_geometry']=='signed_Mobius'
            assert source['closed_source_not_point_sample_extrapolation']
            assert source['nonlinear_original_coefficient_precedes_phase_union']
            assert all(ep(interval(c,b))==(0,1) for b in source['actual_phase_boxes'])
            assert source['actual_phase_left_origin']['positive_original_origin_offset']['strictly_positive']
            for row in source['whole_cell_closed_source_terms']:
                for term in row['terms']:
                    if term['normalized_pressure_late_error_Pstar_power'] is not None:
                        assert term['normalized_pressure_late_error_Pstar_power']==-1
                        assert mp.mp.make_mpf(tuple(term['positive_late_error_finite_budget_upper']['exact_mpf_tuple']))>0
                        late_terms+=1
            for key in totals:
                mass=interval(c,cell['original_positive_own_rate_masses'][key]);assert ep(mass)[0]>0
                totals[key]+=mass
            for record in source['original_inverse_and_Z_piece_records']:
                inverse=record['inverse'];assert ep(interval(c,inverse['coordinate_interval']))==(0,1)
                proof=record['ordinary_Z']['bounded_original_implicit_product']
                assert proof['original_identity']=='t*psi_Z=-T2_Z*t/(1+t^2)'
                assert not proof['selected_endpoint_midpoint_or_source_cap_value']
        assert late_terms>0
        for key,rate in points.exact.RATES.items():
            rr=p.mpf(rate.numerator)/rate.denominator if isinstance(rate,Fraction) else p.mpf(str(rate))
            expected=p.quad(lambda y:p.exp(rr*y),[-5,0])
            contains(totals[key],expected,p.mpf('1e-95'));mass_checks+=1
        contains(totals['p'],c.mpf(5))
        widths.append(dict(cells=count,finite_normalized_integral_magnitude_bounds={
            key:{jet:max(abs(v) for v in ep(saved_value(c,row,owner.atlas.bases)))
                for jet,row in level['actual_finite_N_five_density_contribution_enclosures'][key].items()
                if not any(row['formal_positive_scale']['source_exponents'][:2])}
            for key in ('m','h','k','e','p')}))
        # m/k Z retain L^-1; this original L interval is strictly positive
        # and within1+the positive directed delta budget of unit size.
        for key in ('m','k'):
            row=level['actual_finite_N_five_density_contribution_enclosures'][key]['Z']
            scale=row['formal_positive_scale'];assert scale['source_exponents']==[0,0,-1,0]
            assert max(abs(v) for v in ep(interval(c,scale['additional_log_interval'])))<10
        assert level['incoming_global_correction_histories_not_assumed_zero']
        assert level['original_P0_kept_separate_and_not_reset'] and level['original_pressure_rate_zero_memory_retained']
    return dict(passed=True,independent_positive_own_rate_mass_integral_comparisons=mass_checks,
        exact_zero_rate_pressure_memory_mass=5,pressure_datum_C0_Z_identities=2,
        exact_pressure_affine_template_rows=linear,whole_window_levels=[4,16],
        source_phase_pressure_late_and_inverse_contracts_preserved=True,
        finite_Z_integral_ranges_no_astronomical_implicit_product_artifact=True,
        normalized_ranges=widths)


def run():
    begin=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['signed_implicit_Z_products_bounded_before_multiplication']
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
           'current_whole_N_selected',*current.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[k] is False for k in flags)
    owner=current.OriginalReferenceWholeCells(Z='.37')
    with mp.workdps(owner.ctx.dps+40):
        checks=dict(original_factor_identity_and_ordering=exact_factor_and_ordering(owner),
            signed_implicit_product_and_integral_references=implicit_product_references(owner),
            issued_factor_rebase_and_whole_cell=live_rebase_and_whole_cell(owner),
            full_original_window_rate_and_pressure_contracts=whole_window_contracts(owner,manifest))
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-begin,
        scope='Original fixed nonzero-Z Rh_reference whole-window actual-N C0/Z five-density contributions. Exact source factor correlation and bounded implicit product, finite-unit defining-integral diagnostics, live point/whole-cell and positive own-rate memory checks. Not whole-Z terminal closure, all17/24 controls, selected global-N or recursive field.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original reference whole-cell integrals: factors, implicit product, source and rate memory PASS',flush=True)
    return report


if __name__=='__main__':run()
