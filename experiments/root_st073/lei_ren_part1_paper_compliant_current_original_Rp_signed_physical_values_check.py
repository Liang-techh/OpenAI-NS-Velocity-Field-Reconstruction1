"""Current source-scale identity and independent signed-sum enclosure checks."""
import copy
from dataclasses import replace
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_signed_physical_values as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError):return True
    raise AssertionError('Unproved signed source input admitted')


@source_precision
def run(owner=None,deliveries=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    assert not any(raw[key] for key in current.GATES+current.OPEN)
    if owner is None or not deliveries:raise ValueError('Live candidate and original source deliveries required')
    assert not owner.acceptance_loaded and all(owner.assert_graph().values())
    counts=dict(rows=0,literal_groups=0,canonical_groups=0,numeric_rows=0,relative_width_satisfied_rows=0,
        retained_positive_ratio_tails=0)
    results={};all_proofs=0
    for name,delivery in deliveries.items():
        out=owner.evaluate(delivery);results[name]=out
        assert current.correlated.report(out)==raw['actual_current_signed_physical_values'][name]
        source=owner.before.before.physical_rows(delivery);reader=owner.before.reader(delivery)
        for section,components in out['physical_value_rows'].items():
            for component,rows in components.items():
                original=source[section][component]
                original=original if section=='Cartesian_spatial_rows' else {'dt':original}
                assert set(rows)==set(original)
                for label,value in rows.items():
                    expected=original[label];assert value['physical_derivative']==expected.derivative
                    partition={}
                    for group in expected.groups():
                        log=owner.graph.add(*[ref for _,ref in group['log_scale_parts']])
                        key=(group['scale_units'],group['exact_source_scale_powers'],tuple(sorted(reader.polynomial(log).items())))
                        partition.setdefault(key,[]).append(group)
                    assert len(value['signed_canonical_scale_groups'])==len(partition)
                    for report,((units,powers,polynomial),parts) in zip(value['signed_canonical_scale_groups'],partition.items()):
                        canonical=current.correlated.box.pulse.radius.FunctionRef(owner.graph,report['exact_scale_log_function'])
                        assert tuple(sorted(reader.polynomial(canonical).items()))==polynomial
                        assert tuple(report['same_source_units'])==units
                        assert report['same_source_scale_powers']==[current.rational_record(v) for v in powers]
                        independent=owner.ctx.mpf(0)
                        for group in parts:
                            independent+=group['signed_coefficient'][0]
                            log=owner.graph.add(*[ref for _,ref in group['log_scale_parts']])
                            assert not reader.polynomial(owner.graph.sub(log,canonical))
                        assert report['signed_coefficient_enclosure']._mpi_==independent._mpi_
                        for proof_node in report['exact_scale_identity_proofs']:
                            proof=owner.graph.nodes[proof_node]
                            assert proof['operation']=='exact_original_Rp_canonical_physical_positive_scale_identity'
                            assert tuple(proof['identical_source_units'])==units
                            assert proof['identical_source_scale_powers']==[current.rational_record(v) for v in powers]
                            assert proof['amplitude_identity']==owner.before.proof.node
                            assert proof['same_physical_component']==component
                            a=current.correlated.box.pulse.radius.FunctionRef(owner.graph,proof['original_scale_log_function'])
                            b=current.correlated.box.pulse.radius.FunctionRef(owner.graph,proof['canonical_scale_log_function'])
                            assert not reader.polynomial(owner.graph.sub(a,b));all_proofs+=1
                    for ratio in value['ratio_terms']:
                        reference=current.correlated.box.pulse.radius.FunctionRef(owner.graph,value['exact_reference_log_scale_function'])
                        group=current.correlated.box.pulse.radius.FunctionRef(owner.graph,ratio['source_log_scale'])
                        difference=owner.graph.sub(group,reference)
                        assert difference.node==ratio['exact_log_ratio_function']
                        assert reader.at(difference)._mpi_==ratio['directed_log_ratio']._mpi_
                        method=ratio['method']['kind']
                        if method=='exact_common_scale':assert not reader.polynomial(difference)
                        elif method=='retained_positive_tail':
                            assert current.ends(reader.at(difference))[1]<=-current.EXP_LOG_LIMIT
                            assert current.ends(ratio['directed_ratio_enclosure'])[1]>0
                            assert ratio['method']['tail_not_replaced_by_zero']
                    if value['ordinary_numeric_materialized']:
                        assert value['exact_zero_enclosure'],'Original-scale nonzero expansion was not observed in these source cases'
                        assert current.ends(value['ordinary_numeric_enclosure'])==(0,0)
                    if component=='p':assert value['current_P0_numeric_inclusion_pending']
                    assert value['actual_source_and_original_physical_operators_retained']
        for key in counts:counts[key]+=out['delivery_counts'][key]
        assert out['original_astronomical_scale_not_replaced'] and not any(out[key] for key in current.GATES+current.OPEN)
    # Numeric fixtures test the summation algorithm, not source reconstruction.
    c=owner.ctx;g=owner.graph;reader=owner.before.reader();ind=MPIntervalContext();ind.dps=400
    fixtures={}
    for label,terms in {'different_scales':[(3,'9/8'),(-2,'-3/8')],
                        'negative_sum':[(1,'-2'),(0,'1/5')],
                        'equal_scale_exact_cancel':[(7,'2'),(7,'-2')]}.items():
        groups=[];true=ind.mpf(0);symbolic=s.Integer(0)
        for scale,coefficient in terms:
            ref=g.constant(scale);q=Fraction(coefficient)
            box=c.mpf(q.numerator)/q.denominator
            groups.append(dict(log_function=ref,log_bound=reader.at(ref),coefficient=box))
            true+=(ind.mpf(q.numerator)/q.denominator)*ind.exp(ind.mpf(scale))
            symbolic+=s.Rational(q.numerator,q.denominator)*s.exp(scale)
        # An independent interval sum loses the exact same-scale correlation.
        # Prove the zero algebraically before demanding a zero enclosure.
        if symbolic==0:true=ind.mpf(0)
        result=current.enclose_factored_sum(g,c,reader,groups,Fraction(1,10**8))
        assert result['ordinary_numeric_materialized']
        assert current.amplitude.contains(result['ordinary_numeric_enclosure'],true)
        assert result['requested_relative_width_satisfied'];fixtures[label]=result
    huge=owner.before.radius.logRp;bound=reader.at(huge)
    result=current.enclose_factored_sum(g,c,reader,[dict(log_function=huge,log_bound=bound,coefficient=c.mpf(2)),
        dict(log_function=huge,log_bound=bound,coefficient=c.mpf(-2))],Fraction(1,10**8))
    assert result['exact_cancellation_before_scale_expansion'] and current.ends(result['ordinary_numeric_enclosure'])==(0,0)
    fixtures['exact_cancel_at_actual_astronomical_source_scale']=result
    uncertain=current.signed_log_enclosure(c,c.mpf([-1,1]),c.mpf(3))
    assert uncertain['sign'] is None and uncertain['log_magnitude_lower']=='-infinity'
    tail,method=current.ratio_enclosure(c,c.mpf([-1200,-1100]))
    assert method['tail_not_replaced_by_zero'] and current.ends(tail)[1]>0
    assert current.amplitude.contains(tail,ind.exp(-1100))
    unavailable,method=current.ratio_enclosure(c,c.mpf([-2000,2000]))
    assert unavailable is None and method['kind']=='unresolved_ratio'
    # Unit/power separation regression. These edited metadata terms are only
    # a private schema fixture; the public caller accepts live sources only.
    sample=next(iter(source['Cartesian_spatial_rows']['ux'].values())).terms[0]
    sample=replace(sample,log_scale_parts=(('finite_fixture',g.constant(3)),))
    other=replace(sample,source_row=replace(sample.source_row,powers=sample.source_row.powers+(Fraction(1),)))
    row=replace(next(iter(source['Cartesian_spatial_rows']['ux'].values())),terms=(sample,other))
    ids=dict(original_physical_inverse_identity=out['original_physical_inverse_identity'],
        original_native_inverse_identity=out['original_native_inverse_identity'])
    assert len(owner._canonical_groups(row,reader,ids))==2
    component_views={}
    for case_name,delivery in deliveries.items():
        components=owner.velocity_pressure(delivery);component_views[case_name]=components
        assert set(components['values'])=={'u','v','w','p'}
        _,_,request=owner.before.before._validate_delivery(delivery)
        assert components['physical_coordinates']=={key:ref.node for key,ref in request.forward_coordinates.items()}
        for alias,component in (('u','ux'),('v','uy'),('w','uz'),('p','p')):
            assert current.correlated.report(components['values'][alias])==current.correlated.report(
                results[case_name]['physical_value_rows']['Cartesian_spatial_rows'][component]['x0_y0_z0'])
        assert not components['unrestricted_physical_point_API'] and not components['full_certified_physical_accuracy']
    bad=dict(copied_source_delivery=rejected(lambda:owner.evaluate(copy.copy(next(iter(deliveries.values()))))),
        caller_value_oracle=rejected(lambda:owner.evaluate({'u':[0,0]})),
        invalid_accuracy=rejected(lambda:owner.evaluate(next(iter(deliveries.values())),'0')),
        caller_accuracy_interval=rejected(lambda:owner.evaluate(next(iter(deliveries.values())),c.mpf(['.0001','.001']))))
    assert all(bad.values()) and all(owner.assert_graph().values())
    receipt=dict(all_passed=True,source_family=owner.family_record,actual_source_counts=counts,
        exact_actual_source_positive_scale_identity_proofs=all_proofs,
        actual_signed_common_scale_values=results,independent_400_digit_finite_sum_fixtures=fixtures,
        actual_supported_family_uvw_pressure_component_interfaces=component_views,
        exact_cancellation_at_original_astronomical_scale_checked=True,
        sign_crossing_keeps_unbounded_lower_log=True,strict_positive_ratio_tail_retained=True,
        nonmatching_units_or_powers_not_merged=True,invalid_inputs_rejected=bad,
        original_nonzero_values_remain_factored_not_ordinary_numbers=True,
        P0_numeric_source_inclusion_and_full_physical_accuracy_remain_open=True,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt),indent=2)+'\n',
        encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_SIGNED_PHYSICAL_VALUES actual signed sums and positive tails',flush=True)
    return receipt


if __name__=='__main__':run()
