"""Independent Decimal arithmetic, original graph binding and phase errors.

Only the new phase theorem is checked. No ancestor quadratures are rerun,
and no interval cap is read as a field point or selected global frequency.
"""
from decimal import Decimal,localcontext
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_all_chart_point_phase as current
import lei_ren_part1_paper_compliant_current_native_Rc_function_transport_check as transport_check
import lei_ren_part1_paper_compliant_current_original_O2_radius_phase_points as O2


def interval(c,record):
    p=mp.mp.clone();p.dps=c.dps+40
    return c.mpf((p.make_mpf(tuple(record['lower_exact_mpf_tuple'])),p.make_mpf(tuple(record['upper_exact_mpf_tuple']))))


def exact_selected_fraction(bits,multiplier):
    sign,man,exponent,_=bits
    # Independent modular numerator/denominator arithmetic. A positive large
    # exponent only enters modular exponentiation, never integer expansion.
    num=(-1 if sign else 1)*man*multiplier.numerator;den=multiplier.denominator
    if exponent < -8192:raise ValueError('Microscopic binary input must stay factored')
    if exponent<0:den*=2**(-exponent);rem=num%den
    else:rem=(num%den)*pow(2,exponent,den)%den
    return Fraction(rem,den)


def decimal_references(owner,manifest):
    c=MPIntervalContext();c.dps=850;p=mp.mp.clone();p.dps=850
    samples=[*manifest['original_all17_true_point_phase_records'].values(),
        *manifest['large_explicit_N_point_phase_records'].values(),
        *manifest['original_source_expression_endpoint_records'].values()]
    tiny_components=0;selected_components=0;max_width=p.mpf(0)
    sc=c.mpf(owner.binder.sc)
    assert current.ep(sc)[0]==current.ep(sc)[1] and current.ep(sc)[0]>0
    with localcontext() as dc:
        dc.prec=800;one=Decimal(1)
        P=Decimal(40).exp()+11;W=60*Decimal(1000).ln()+240*P
        L100=Decimal(100).ln()-Decimal(4).ln()+1000
        L110=Decimal(110).ln()-Decimal(4).ln()+1000
        for record in samples:
            chart=record['chart'];N=record['explicit_candidate_N'];cr=record['coordinate'];q=Fraction(cr['exact_coefficient'])
            dq=Decimal(q.numerator)/Decimal(q.denominator)
            iq=c.mpf(q.numerator)/q.denominator
            if cr['kind']=='selected_sc_multiple':x=None;ix=sc*iq
            elif cr['kind']=='original_power_offset':x=dq/W;ix=iq/owner.analytic(owner.symbols['Tw'],c)
            elif cr['kind']=='original_patch_log_offset':x=dq.exp();ix=c.exp(iq)
            else:x=dq;ix=iq
            micro={'hbB':-sc/2};dyadic={};regular=14*P+L110
            if chart in ('bridge_first','bridge_second'):
                regular=Decimal(0)
                micro={'hbB':sc*(iq-c.mpf('.5')) if cr['kind']=='selected_sc_multiple' else ix-sc/2}
            elif chart=='bridge_macro':regular=(4*P+L100)*x;micro={'hbB':2*(1-ix)-sc/2}
            elif chart in ('switch_first','switch_second'):regular=4*P+L100;micro['hbS']=ix
            elif chart=='switch_power':regular=4*P+L100+(Decimal(110).ln()-Decimal(100).ln())*x;micro['hbS']=2*(1-ix)
            elif chart=='reshape':regular=4*P+L110;dyadic={'T':q}
            elif chart=='inner_reference':regular=(4+10*x)*P+L110-8*x;dyadic={'logC':10*q,'T':1-q}
            else:
                dyadic={'logC':Fraction(10)}
                if chart=='axial_restore':regular+=x-8
                elif chart=='restore_buffer':regular+=x
                elif chart=='actual_patch':regular+=(-6+dq if cr['kind']=='original_patch_log_offset' else -6+x.ln())
                elif chart in ('Rh_reference','O2_slope'):regular+=x
                elif chart=='O2_axial':regular+=(40*x).exp()
                elif chart=='O2_buffer':regular+=Decimal(40).exp()+x
                elif chart=='O3_slope_mu':regular+=P+x
                elif chart=='O3_power':regular+=P+1+(dq if cr['kind']=='original_power_offset' else W*x)
                else:raise AssertionError('Missing independent original radius chart')
            reference=regular*N
            for name,coefficient in dyadic.items():
                bits=current.ep(owner.binder.fixed[name])[0]._mpf_;fraction=exact_selected_fraction(bits,coefficient*N)
                reference+=Decimal(fraction.numerator)/Decimal(fraction.denominator)
            reference-=reference//one;ref=p.mpf(str(reference))
            boxes=[interval(c,box) for box in record['actual_phase_directed_boxes']]
            # This is the exact regular/source-selected part with zero offset.
            # Directed signed microscopic budgets include zero and the actual
            # source offset; the original widths are never evaluated as zero.
            assert any(current.ep(box)[0]<=ref<=current.ep(box)[1] for box in boxes),(chart,N)
            width=sum((p.mpf(current.ep(box)[1])-p.mpf(current.ep(box)[0]) for box in boxes),p.mpf(0))
            assert width<=p.mpf(10)**(-record['requested_phase_digits'])
            max_width=max(max_width,width)
            assert not record['periodic_projection']['full_period']
            for item in record['original_regular_component_proofs']:
                if item['method']!='exact_selected_dyadic_modulus':continue
                coefficient=Fraction(item['exact_multiplier']['numerator'],item['exact_multiplier']['denominator'])
                assert coefficient==dyadic[item['component']]*N
                exact=exact_selected_fraction(item['source_exact_mpf_tuple'],coefficient)
                assert item['source_exact_mpf_tuple']==list(current.ep(owner.binder.fixed[item['component']])[0]._mpf_)
                assert item['exact_fraction']==dict(numerator=exact.numerator,denominator=exact.denominator)
                assert not item['huge_integer_materialized'];selected_components+=1
            seen=set()
            for item in record['original_signed_microscopic_phase_components']:
                name='hbB' if item['original_width_kind']=='bridge' else 'hbS';seen.add(name)
                coefficient=micro[name];cc=interval(c,item['coefficient_enclosure']);lo,hi=current.ep(cc)
                clo,chi=current.ep(coefficient)
                assert clo<=hi and lo<=chi
                budget=interval(c,item['numerical_phase_error_enclosure']);blo,bhi=current.ep(budget)
                assert blo<=0<=bhi and item['actual_width_not_zeroed'] and item['actual_width_not_materialized']
                if (clo,chi)!=(0,0):
                    assert not item['exact_zero_coefficient']
                    logwidth=owner.binder.loghB if name=='hbB' else owner.binder.loghS
                    assert current.ep(interval(c,item['original_positive_width_log']))==current.ep(c.mpf(logwidth))
                    independent=c.ln(N)+c.mpf(logwidth)+c.ln(c.mpf(max(abs(lo),abs(hi))))
                    bound=interval(c,item['original_N_width_coefficient_phase_log_upper'])
                    assert current.ep(independent)[0]<=current.ep(bound)[1] and current.ep(bound)[0]<=current.ep(independent)[1]
                    assert current.ep(bound)[1]<current.ep(c.ln(max(abs(blo),abs(bhi))))[0]
                    if clo>=0:assert blo==0 and bhi>0
                    elif chi<=0:assert blo<0 and bhi==0
                    else:assert blo<0<bhi
                tiny_components+=1
            assert seen=={name for name,value in micro.items() if current.ep(value)!=(0,0)}
            assert record['source_family']==owner.family and record['source_graph_sha256']==owner.source_graph_sha256
            assert record['original_width_Jacobians_unchanged'] and record['original_radius_Z_derivative_exact_zero']
            assert record['source_phase_not_an_independent_free_angle'] and record['source_field_values_not_selected_from_caps']
            assert not record['full_original_source_point_or_integral_oracle_installed']
    assert len(samples)==24 and len(manifest['large_explicit_N_point_phase_records'])==3
    assert all(row['adaptive_interval_digits']>260 and row['explicit_candidate_N'].bit_length()==1201
        for row in manifest['large_explicit_N_point_phase_records'].values())
    zero=manifest['original_source_expression_endpoint_records']['zero_inlet']
    assert zero['exact_original_radius_minus_inlet_expression']==sy.srepr(sy.Integer(0))
    assert not zero['original_signed_microscopic_phase_components']
    assert sy.diff(owner.maps['bridge_first'],owner.x)==owner.symbols['hbB']
    return dict(passed=True,independent_Decimal_digits=800,independent_point_phase_queries=len(samples),
        original_selected_dyadic_component_checks=selected_components,signed_retained_source_width_error_checks=tiny_components,
        maximum_total_periodic_cover_width=max_width,large_1201_bit_N_adaptive_queries=3,
        original_selected_sc_remains_exact_binary_without_rational_expansion=True,
        exact_zero_inlet_offset_keeps_original_nonzero_width_Jacobian=True)


def graph_bindings(owner):
    manifest=json.loads((current.HERE/current.functions.NAME).read_bytes())
    g=current.functions.FunctionTransportGraph()
    for index,row in enumerate(manifest['function_graph_nodes']):
        assert g.node(**row).node==index
    N=next(i for i,row in enumerate(g.nodes) if row['operation']=='shared_positive_integer')
    parameters={row['name']:current.functions.FunctionRef(g,i) for i,row in enumerate(g.nodes)
        if row['operation']=='original_source_parameter'}
    built=dict(graph=g,N=current.functions.FunctionRef(g,N),parameters=parameters,
        source_family=owner.family,source_graph_sha256=owner.source_graph_sha256)
    count=0;before=len(g.nodes)
    saved=json.loads((current.HERE/current.NAME).read_bytes())['original_all17_true_point_phase_records']
    for row in g.nodes:
        if row['operation']!='original_function_graph':continue
        coordinate={'original_power_offset':'2'} if row['chart']=='O3_power' else saved[row['chart']]['coordinate']['exact_coefficient']
        value=owner.phase_for_source(row,built,coordinate=coordinate,N=257)
        assert value['chart']==row['chart'];count+=1
    assert count>0 and len(g.nodes)==before
    selected=next(row for row in g.nodes if row['operation']=='original_function_graph')
    endpoint=next(row for row in g.nodes if row['operation']=='original_function_graph' and row['chart']=='O3_power')
    rejected=0
    bad_calls=[lambda:owner.phase_for_source(dict(selected),built,coordinate='.5',N=257),
        lambda:owner.phase_for_source(selected,{**built,'source_graph_sha256':'wrong'},coordinate='.5',N=257),
        lambda:owner.phase_for_source(selected,built,coordinate='.5',N=159),
        lambda:owner.phase_for_source(selected,built,coordinate='.5',N=1),
        lambda:owner.phase_for_source(endpoint,built,coordinate={'original_power_offset':'1'},N=257),
        lambda:owner.evaluate(chart='invalid',coordinate='.5',N=257),
        lambda:owner.evaluate(chart='O2_slope',coordinate='-1',N=257),
        lambda:owner.evaluate(chart='O2_slope',coordinate=capped_interval(),N=257),
        lambda:owner.evaluate(chart='O2_slope',coordinate='.5',N=True),
        lambda:owner.evaluate(chart='O2_slope',coordinate='.5',N=0),
        lambda:owner.evaluate(chart='O2_slope',coordinate='.5',N=1<<4096),
        lambda:owner.evaluate(chart='O2_slope',coordinate='.5',N=257,decimal_digits=39)]
    for call in bad_calls:
        try:call()
        except (ValueError,TypeError):rejected+=1
    assert rejected==len(bad_calls) and len(g.nodes)==before
    for diagnostic_N in (1,159):
        diagnostic=owner.evaluate(chart='O2_slope',coordinate='.537',N=diagnostic_N)
        assert diagnostic['candidate_N_not_global_frequency_admission']
        assert not diagnostic['full_original_source_point_or_integral_oracle_installed']
    return dict(passed=True,accepted_original_function_phase_rows_checked=count,
        original_endpoint_coordinate_and_N_bound=True,source_graph_not_mutated=True,
        sub160_phase_only_diagnostics_cannot_enter_original_graph=True,
        invalid_source_domain_coordinate_frequency_and_precision_rejections=rejected)


def capped_interval():
    c=MPIntervalContext();c.dps=80;return c.mpf(('.4','.6'))


def wrap_and_O2(owner):
    c=MPIntervalContext();c.dps=120
    wrapped=current.periodic_add(c,[c.mpf(('.999','.9995'))],c.mpf(('.00025','.001')))
    assert not wrapped['full_period'] and len(wrapped['boxes'])==2
    assert any(current.ep(box)[0]<=mp.mpf('.9995')<=current.ep(box)[1] for box in wrapped['boxes'])
    assert any(current.ep(box)[0]<=mp.mpf('.00025')<=current.ep(box)[1] for box in wrapped['boxes'])
    O2owner=O2.OriginalO2RadiusPhasePoints();comparisons=0
    for chart,y in (('O2_slope','.537'),('Rh_reference','-2.337')):
        first=owner.evaluate(chart=chart,coordinate=y,N=257)
        # O2's accepted radius formula also applies to reference offsets. Its
        # public domain guard excludes negative y, so use an equivalent phase
        # y+integer with integer N and check the exact periodic translation.
        other=O2owner.evaluate(y='.663' if chart=='Rh_reference' else y,N=257)
        oldboxes=other['true_original_phase_directed_boxes']
        def scalar(value):
            return mp.make_mpf(tuple(value._mpf_)) if hasattr(value,'_mpf_') else value
        assert any(max(current.ep(a)[0],scalar(b['lower']))<=min(current.ep(a)[1],scalar(b['upper']))
            for a in first['actual_phase_directed_boxes'] for b in oldboxes)
        comparisons+=1
    return dict(passed=True,circular_Minkowski_union_keeps_wrap_seam=True,
        accepted_O2_true_radius_recipe_overlap_comparisons=comparisons,
        reference_offset_integer_translation_exact=True),O2owner.hashes


def run(owner=None):
    began=time.monotonic();fresh=owner is None
    if owner is None:
        bridge,_=current.phase.native.inlet.native_bridge_owner()
        with current.phase.native.inlet.CheckedSourceRuntime():
            owner=current.OriginalAllChartPointPhase(current.phase.NativeSpatialPhase(current.phase.native.NativeGenericSourcePackets(bridge)))
    report=json.loads((current.HERE/current.NAME).read_bytes())
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('full_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.phase.packets.OPEN)
    assert all(report[key] is False for key in flags)
    geometry=transport_check.geometry_reference()
    references=decimal_references(owner,report);bindings=graph_bindings(owner);overlaps,O2hashes=wrap_and_O2(owner)
    hashes={**report['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
        Path(transport_check.__file__).name:current.sha(Path(transport_check.__file__).name),Path(O2.__file__).name:current.sha(Path(O2.__file__).name)}
    for name,digest in O2hashes.items():
        assert current.sha(name)==digest,name
        if name in hashes:assert hashes[name]==digest,name
        hashes[name]=digest
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        independent_original_radius_and_Jacobian_geometry=geometry,independent_all17_Decimal_and_signed_error_references=references,
        accepted_original_function_source_binding=bindings,periodic_wrap_and_original_O2_overlap=overlaps,
        source_owner_reused_only_when_live=True,standalone_fresh_original_owner_constructed=fresh,
        ancestor_producers_or_quadratures_reexecuted=False,**dict.fromkeys(flags,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='All17 original exact-point phase service with directed errors only. Source leaf C0/Z/integral callbacks, control values, global N and actual recursion remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(result),indent=2).encode()+b'\n')
    print('Original17 chart true point phase: independent Decimal, source graph, width errors and wraps PASS',flush=True)
    return result


if __name__=='__main__':run()
