"""Focused directed-arithmetic and original live-root checks."""
import json
import time
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_signed_input_enclosures as current

packets=current.packets;HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha


def require(value,message):
    if not value:raise ArithmeticError(message)


def contains(outer,inner):
    a,b=current.ep(outer);c,d=current.ep(inner)
    return a<=c and d<=b


def arithmetic_check(c):
    """Independent bounded exact values exercise signed cancellation and caps."""
    bases=tuple(c.mpf(v) for v in ('1e40','0','0','0','0'))
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    value=lambda coefficient,power=0,offset=0:current.ScaledEnclosure(
        current.FormalScale(bases,(power,0,0,0,0),offset),coefficient,ledger)
    numerator,denominator=value(-3,1),value(2,1)
    quotient=numerator.positive_divide(denominator,bases[0]+c.ln(2))
    require(contains(quotient.finite_interval(),c.mpf('-1.5')),'Signed huge-factor quotient must cancel exactly')
    require((value(7,1)-value(7,1)).zero,'Same formal source must cancel before exponentiation')
    guarded=False
    try:numerator.finite_interval()
    except ArithmeticError:guarded=True
    require(guarded,'Huge source exponential must remain forbidden')
    uncertain=value((-1,2))
    reciprocal=value(1).positive_divide(uncertain,c.ln(c.mpf('.25')))
    require(contains(reciprocal.finite_interval(),c.mpf(('.5','4'))),'Positive theorem must enclose every admissible reciprocal')
    # An omitted negative tail must not become exact zero or a positive-only cap.
    tail=value(1)+value(-3,offset='-100000')
    require(current.ep(tail.coefficient)[0]<1 and current.ep(tail.coefficient)[1]>=1,
        'Directed negative tail must be retained')
    # Independent wide log intervals exercise coordinate rescaling. The true
    # values at both endpoints remain in its returned signed interval.
    mixed=value(2,offset=c.mpf((0,10000)))+value(-1,offset=c.mpf((0,9999)))
    require(mixed.ledger['directed_independent_log_rescalings']>0,'Wide log rescaling path must execute')
    require(current.ep(mixed.coefficient)[0]<0 and current.ep(mixed.coefficient)[1]>0,
        'Rescaled wide signed box must retain both signs')
    return dict(passed=True,signed_huge_factor_cancellation=True,
        huge_exponential_materialization_rejected=True,
        positive_crossing_denominator_intersection_checked=True,
        negative_exponential_tail_not_discarded=True,independent_wide_log_rescaling_checked=True,
        arithmetic_ledger=ledger)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['live_native_chart_count']==17,'All17 live signed source records required')
    require(not any(saved.get(k) for k in packets.OPEN),'Root-stage global gates must remain false')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed native root input: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))
        fixture=arithmetic_check(owner.ctx)
        count=0;regions={};amplitude_checks=0;O3_checks=0
        for chart,record in saved['current_native_signed_source_root_records'].items():
            provenance=record['source_provenance']
            Z=packets.interval(owner.ctx,provenance['Z_box']);v=packets.interval(owner.ctx,provenance['coordinate_box'])
            live=owner.query(chart,Z,v);fresh=packets.encode(live['record'])
            require(fresh==record,'Stored signed roots differ from the current live source: '+chart)
            roots=live['roots'];positive=owner.decode(owner.inventory[chart]['actual_positive_denominator_theorem'])
            for name,row in roots.items():
                require(len(row)==6,'All signed y2/Z1 root orders required')
                for value in row.values():
                    require(not value.record()['point_value_selected'],'Interval endpoint cannot define a field value')
                    count+=1
            for name,key in (('E','log_E_positive_lower'),('a','log_actual_a_positive_lower')):
                lower,upper=current.ep(roots[name]['y0_Z0'].coefficient)
                require(lower>0,'Actual positive source root must retain a positive enclosure')
            if provenance.get('amplitude_adapter'):
                proof=provenance['amplitude_adapter']
                require(proof['positive_source_cap_not_called'] and proof['original_row_formula_unchanged'],
                    'Actual frozen amplitude must remain a defining source factor')
                require(set(proof['positive_exponential_calls'])=={1,2},'Original amplitude and square required')
                amplitude_checks+=1
            if chart in ('inner_reference','axial_restore','restore_buffer'):
                # Original logU_y=1/10 throughout these three charts.
                require(contains(roots['a']['y0_Z0'].finite_interval(),owner.ctx.mpf('.8')),
                    'Original restore shear a=4/5 must be enclosed')
            if chart.startswith('O3_'):
                mu=packets.interval(owner.ctx,positive['actual_positive_mu'])
                sigma=current.sigma_jets(owner.ctx,v) if chart=='O3_slope_mu' else None
                for j,k in current.signed.ORDERS:
                    expected=0 if k else (2*mu*sigma[j]*current.math.factorial(j) if sigma is not None else 2*mu if j==0 else 0)
                    actual=roots['kappa_minus2']['y%d_Z%d'%(j,k)]
                    expected_value=current.ScaledEnclosure(current.FormalScale(actual.scale.bases),expected,actual.ledger)
                    require(packets.encode(actual.record())==packets.encode(expected_value.record()),
                        'Every O3 ordinary excess derivative must use the unrounded source')
                    O3_checks+=1
            require(provenance['formal_radius_correlations_used_in_arithmetic'] is False,
                'Conservative radius covers must not be promoted to exact radius correlation')
            regions[chart]=dict(root_count=sum(len(row) for row in roots.values()),
                branch=live['record']['branch_against_actual_eta'],arithmetic_ledger=live['record']['numerical_arithmetic_ledger'])
            print('Live signed enclosure checked:',chart,flush=True)
    require(count==714 and amplitude_checks==5 and O3_checks==12,'Expected original root/amplitude/O3 coverage')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        live_signed_root_enclosures_checked=count,live_native_charts_checked=len(regions),
        original_formal_amplitude_adapters_checked=amplitude_checks,unrounded_O3_excess_ordinary_derivatives_checked=O3_checks,
        focused_directed_arithmetic_check=fixture,regions=regions,
        source_caps_or_midpoints_used_as_field_values=False,
        actual_phase_inverse_or_changed_integral_admission=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Live signed root enclosure backend PASS:',count,'roots across17 charts',flush=True)
    return result


if __name__=='__main__':run()
