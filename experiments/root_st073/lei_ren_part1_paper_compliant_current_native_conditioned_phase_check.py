"""Independent scalar references, narrow-peak cases and fresh native replay."""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_conditioned_phase as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def contains(a,b):
    lo,hi=ep(a);x,y=ep(b);return lo<=x and y<=hi


def fixture(c,a,b,p2,q,E='1.3'):
    bases=tuple(c.mpf(0) for _ in range(5));ledger=dict(directed_small_exponential_tails=0,
        positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    scalar=lambda v:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),v,ledger)
    a,b,p2,q,E=[scalar(v) for v in (a,b,p2,q,E)]
    t0=(-b).positive_divide(a,c.ln(a.coefficient))
    return dict(q=q,roots={k:{(0,0):v} for k,v in dict(a=a,b=b,p2=p2,t0=t0,E=E).items()}),scalar


def scalar_reference_checks(c):
    scal=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    comparisons=0;crosssign=0;small=0;flat=0
    for a,b,p1,p2 in (('.8','.2','5','.2'),('.8','-.2','5','-.2'),('.8','.2','5','0'),('2.019','0','5','0')):
        ref=original.GenericShearLoop(scal,a=a,b=b,p1=p1,p2=p2,Utheta='1.3')
        source,_=fixture(c,a,b,p2,mp.nstr(ref.q,105))
        loop=current.ConditionedPhase(source,c.ln(c.mpf(mp.nstr(scal.d_star,105))))
        crosssign+=loop.geometry=='signed_Mobius';small+=loop.geometry=='small_r_series';flat+=loop.flat
        for phase in ('.137','.337','.663','.863'):
            got=loop.evaluate(phase)
            require(got['status']=='enclosed','Bounded source fixture inverse absent')
            selected=got['selected_inverse'];box=selected['coordinate_interval'];low,high=ep(box)
            # Compare the directed inverse bracket against an independent
            # original scalar evaluator at both endpoints, with its roundoff
            # explicitly allowed. No scalar midpoint is used as a field.
            values=[]
            for x in (low,high):
                fraction=ref.ctx.mpf(mp.nstr(x,110))
                if selected['chart']=='E':
                    if fraction in (0,ref.ctx.mpf('.5'),1):psi=2*ref.ctx.pi*fraction
                    else:
                        rho=ref.one_minus_abs_r
                        plus,minus=(rho,2-rho) if ref.r>0 else (2-rho,rho)
                        psi=2*ref.ctx.atan2(plus*ref.ctx.sin(ref.ctx.pi*fraction),minus*ref.ctx.cos(ref.ctx.pi*fraction))
                        if psi<0:psi+=2*ref.ctx.pi
                else:psi=2*ref.ctx.pi*fraction
                values.append(ref.phase_at_angle(psi))
            target=ref.ctx.mpf(phase);error=ref.ctx.mpf('1e-90')
            require(values[0]<=target+error and values[1]>=target-error,'Original scalar phase root not inside directed bracket')
            reference=ref.evaluate(phase)
            for key,rkey in (('A','A'),('B_over_Pstar','B')):
                value=loop.primitives(box,selected['chart'])[key]
                covered=current.bounded_value(value)+c.mpf(('-1e-85','1e-85'))
                require(contains(covered,c.mpf(mp.nstr(reference[rkey],105))),'Original scalar primitive not enclosed: '+key)
            comparisons+=1
    require(crosssign==2 and small==1 and flat==1,'Expected signed/small/flat reference coverage')
    return dict(passed=True,independent_original_scalar_phase_and_AB_comparisons=comparisons,
        opposite_signed_r_checked=True,exact_p2_zero_checked=True,flat_checked=True,
        reference_roundoff_allowed='1e-85 primitive /1e-90 phase',fixtures_are_not_native_source_values=True)


def extreme_checks(c):
    counts=0
    for sign in (-1,1):
        source,scalar=fixture(c,'.8','0',sign,c.sqrt(c.mpf('.75')))
        # The gigantic logarithm cannot be exponentiated. p2's formal factor
        # controls u while q/t0 remain ordinary, well-conditioned inputs.
        source['roots']['p2'][(0,0)]=current.prior.ScaledEnclosure(
            current.prior.FormalScale(source['q'].scale.bases,offset='1e40'),sign,source['q'].ledger)
        loop=current.ConditionedPhase(source,c.ln(c.mpf('.005')))
        require(loop.geometry=='signed_Mobius' and ep(loop.rho.coefficient)[0]>0,'Huge signed u must retain positive factored rho')
        for phase in ('.137','.337','.663','.863'):
            got=loop.evaluate(phase)
            require(got['status']=='enclosed','Huge-log phase must execute')
            image=got['selected_inverse']['phase_image']
            require(contains(image,c.mpf(phase)) and ep(image)[1]-ep(image)[0]<mp.mpf('1e-9'),'Two-angle inverse failed to resolve huge-log peak')
            counts+=1
    source,scalar=fixture(c,2,0,(-2,3),1)
    source['q']=current.prior.ScaledEnclosure(current.prior.FormalScale(source['q'].scale.bases,offset='-1e40'),1,source['q'].ledger)
    loop=current.ConditionedPhase(source,c.ln(c.mpf('.005')));got=loop.evaluate('.137')
    require(loop.geometry=='small_r_series','Positive tiny q with signed p2 must enter small-r formula')
    values=loop.primitives(got['selected_inverse']['coordinate_interval'])
    require(not values['A'].zero and not values['B_over_Pstar'].zero,'Tiny nonzero source primitives must remain factored')
    # Wide signed u cannot be arbitrarily assigned a sign.
    source,scalar=fixture(c,'.8','0',(-1,1),1)
    loop=current.ConditionedPhase(source,c.ln(c.mpf('.005')))
    require(loop.evaluate('.137')['status']=='requires_signed_source_refinement','Wide signed source box must remain unresolved')
    return dict(passed=True,unmaterializable_signed_u_inverse_queries=counts,
        positive_factored_rho_and_s_retained=True,positive_tiny_q_primitives_not_zeroed=True,
        wide_signed_u_not_given_a_selected_sign=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_source_box_count']==5,'Five actual source boxes required')
    require(not any(saved.get(k) for k in packets.OPEN),'C0 phase backend cannot complete global stages')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed conditioned phase prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeConditionedPhase(current.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))
        fixtures=scalar_reference_checks(owner.ctx);extremes=extreme_checks(owner.ctx);regions={};count=0
        for chart,record in saved['actual_C0_phase_inverse_and_primitive_records'].items():
            p=record['source']['source_provenance'];Z=packets.interval(owner.ctx,p['Z_box']);coordinate=packets.interval(owner.ctx,p['coordinate_box'])
            source,loop=owner.query(chart,Z,coordinate)
            require(packets.encode(source['record'])==record['source'],'Changed live source box: '+chart)
            require(packets.encode(loop.geometry_record())==record['conditioned_geometry'],'Changed live conditioned geometry: '+chart)
            widths=[]
            for phase,previous in record['candidate_phase_queries'].items():
                live=loop.evaluate(phase)
                require(packets.encode(live)==previous,'Changed live phase result: '+chart+' '+phase)
                require(live['status']=='enclosed' and live['inverse_installed_on_this_box'],'Native inverse not available on admitted box')
                image=live['selected_inverse']['phase_image'];require(contains(image,owner.ctx.mpf(phase)),'Inverse phase image lost target')
                width=ep(image)[1]-ep(image)[0];require(width<mp.mpf('1e-9'),'Native phase bracket too wide')
                widths.append(width);count+=1
                if phase in ('0','.5','1') or loop.flat:
                    require(all(v['exact_zero'] for v in live['primitives'].values()),'Endpoint/half-period/flat primitives must be exact zero')
                require(live['free_phase_parameter_not_spatial_phase'] and not live['original_common_N_and_radius_phase_bound'],'Candidate phase cannot become an admitted spatial phase')
            if chart=='O2_buffer':
                require(not record['candidate_phase_queries']['.137']['primitives']['A']['exact_zero'],'Tiny buffer A was incorrectly zeroed')
                require(not record['candidate_phase_queries']['.137']['primitives']['B_over_Pstar']['exact_zero'],'Tiny buffer B was incorrectly zeroed')
            regions[chart]=dict(geometry=loop.geometry,phase_queries=len(widths),max_phase_width=owner.ctx.mpf(max(widths)))
            print('Live conditioned inverse/AB checked:',chart,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        original_native_source_boxes_checked=len(regions),native_phase_inverse_and_AB_queries_checked=count,
        independent_original_scalar_reference_checks=fixtures,focused_extreme_log_checks=extremes,regions=regions,
        spatial_phase_binding_or_slow_derivative_or_integral_admission=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Live conditioned phase inverse/AB PASS:',len(regions),'source boxes;',count,'queries',flush=True)
    return result


if __name__=='__main__':run()
