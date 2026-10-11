"""Independent bounds for actual multi-time point values and fitted changes."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rp_time_observations as current
import lei_ren_part1_paper_compliant_current_original_Rp_remainder_ratios_check as relative_check
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

bg=current.bg
contains=bg.axial.differential.product.contains
rejected=relative_check.stress_check.rejected


@source_precision
def run(owner,value):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not raw[current.GATE] and not any(raw[key] for key in current.OPEN)
    assert raw['source_family']==owner.family_record
    assert raw['input_hashes']==owner.hashes
    assert raw['input_hashes'][Path(current.__file__).name]==current.sha(Path(current.__file__).name)
    view=owner.report(value)
    assert raw['actual_multitime_original_observations']==bg.correlated.report(view)
    entry=owner._fields[id(value)];c=MPIntervalContext();c.dps=800
    deliveries=[record['delivery'] for record,_ in entry[2].values()]
    reader=owner._merged_reader(deliveries,c)
    coefficients={};products_checked=bounds_checked=cone_bounds_checked=0
    for name,sample in entry[1].items():
        record,request=owner._validate_sample(sample);delivery=record['delivery']
        local_reader=relative_check.original_reader(owner,delivery,c)
        assert request.Z!=0 and view['actual_time_samples'][name]['strict_point_cone_certified']
        coordinates=request.forward_coordinates
        identity=owner.graph.add(owner.graph.mul(owner.graph.constant(2),coordinates['loglambda']),
            owner.graph.neg(coordinates['log_tau']),owner.graph.unary('log',owner.graph.sub(owner.graph.one,
                owner.graph.mul(coordinates['Z'],coordinates['Z']))))
        assert not local_reader.polynomial(identity)
        coefficients[name]={}
        source=owner.product.before.source(delivery)
        for label,item in entry[3][name].items():
            if label in ('u','v','w'):
                definitions=dict(u=[(-1,'Ur',bg.CS),(-1-bg.DELTA,'Utheta',-bg.SN)],
                    v=[(-1,'Ur',bg.SN),(-1-bg.DELTA,'Utheta',bg.CS)],w=[(-1-bg.DELTA,'Uz',bg.s.Integer(1))])
                actual=[]
                for term in item['row']['actual_source_product_ledger']:
                    atoms=term['actual_shared_source_product']
                    assert len(atoms)==1 and atoms[0]['ordinary_derivative']==(0,0) and atoms[0]['power']==1
                    assert term['radial_power']==bg.signed.rational_record(Fraction(0))
                    actual.append((bg.s.sympify(term['lambda_exponent']),atoms[0]['source'],
                        bg.s.sympify(term['exact_operator_expression'])))
                assert len(actual)==len(definitions[label])
                for beta,primitive,operator in definitions[label]:
                    matches=[a for a in actual if a[1]==primitive and bg.s.simplify(a[0]-beta)==0]
                    assert len(matches)==1 and bg.s.simplify(matches[0][2]-operator)==0
                assert item['unit']=='velocity'
                if label in ('u','v'):assert item['expected_tau_exponent'] is None
            if item['row'] is not None:
                coefficient,count=relative_check.independent_row(owner,delivery,item['row'],c,
                    source=source,request=request,reader=local_reader)
                assert coefficient is not None;products_checked+=count
            elif label=='omega_z':
                assert item['unit']==owner.arithmetic._require(sample.axial)['unit']=='velocity/length**1'
                theta=source['log_radius_mixed_rows']['Utheta']['y0_Z0']
                coefficient=-c.sqrt(2)*c.mpf(theta.coefficients[0])
                log=owner.graph.add(*(ref for _,ref in theta.log_scale_parts),owner.product.logmu,
                    owner.graph.mul(owner.graph.constant(Fraction(-1,2)),coordinates['logR']),
                    owner.graph.mul(owner.graph.sub(owner.graph.constant(-2),owner.product.delta),coordinates['loglambda']))
                assert not local_reader.polynomial(owner.graph.sub(log,item['scale']))
            else:
                assert label in ('tau','lambda','radius','abs_axial_coordinate','axial_to_radial_coordinate_ratio')
                coefficient=c.mpf(1)
                expected_scale=dict(tau=coordinates['log_tau'],lambda_=coordinates['loglambda'],radius=coordinates['log_r'],
                    abs_axial_coordinate=request.physical_log_abs_z,
                    axial_to_radial_coordinate_ratio=owner.graph.sub(request.physical_log_abs_z,coordinates['log_r']))
                assert not local_reader.polynomial(owner.graph.sub(item['scale'],expected_scale['lambda_' if label=='lambda' else label]))
            assert contains(item['coefficient'],coefficient)
            lo,hi=bg.ends(coefficient);assert (lo>0 if item['sign']==1 else hi<0)
            coefficients[name][label]=coefficient;bounds_checked+=1
        margins=owner.cone._fields[id(sample.cone)][4]
        for label,row in margins.items():
            coefficient,count=relative_check.independent_row(owner,delivery,row,c,
                source=source,request=request,reader=local_reader)
            assert coefficient is not None and bg.ends(coefficient)[0]>0
            products_checked+=count;cone_bounds_checked+=1
        for label,(a,b) in current.relative.PAIRS.items():
            a='Er_total' if a=='Er' else a
            numerator,denominator=entry[3][name][a],entry[3][name][b]
            bound=local_reader.at(owner.graph.sub(numerator['scale'],denominator['scale']))+\
                c.ln(abs(coefficients[name][a]/coefficients[name][b]))
            diagnostic=view['actual_time_samples'][name]['finite_remainder_relative_log_magnitudes'][label]
            assert contains(diagnostic['directed_log_absolute_relative_value'],bound)
            lo,hi=bg.ends(bound)
            if diagnostic['magnitude_strictly_below_one']:assert hi<0
            if diagnostic['magnitude_strictly_above_one']:assert lo>0
    observations={};fits_checked=relative_fits_checked=0
    names=list(entry[1])
    def check_fit(fit,left,right,q,ca,cb):
        difference=owner.graph.sub(right['scale'],left['scale'])
        ref=bg.box.pulse.radius.FunctionRef(owner.graph,fit['exact_scale_log_change'])
        assert not reader.polynomial(owner.graph.sub(ref,difference))
        log_bound=reader.at(difference);coefficient=cb/ca
        assert bg.ends(coefficient)[0]>0
        assert contains(fit['directed_scale_log_change'],log_bound)
        assert contains(fit['actual_coefficient_magnitude_ratio'],coefficient)
        magnitude=log_bound+c.ln(coefficient)
        exponent=magnitude/(c.mpf(q.numerator)/q.denominator)
        assert contains(fit['directed_log_magnitude_ratio'],magnitude)
        assert contains(fit['fitted_tau_exponent'],exponent)
        numeric=relative_check.independent_ratio(c,coefficient,log_bound)
        if numeric is None:assert fit['ordinary_magnitude_ratio'] is None
        else:assert contains(fit['ordinary_magnitude_ratio'],numeric)
        expected=left['expected_tau_exponent']
        if expected is None:assert fit['original_expected_tau_exponent'] is None
        else:
            model=reader.at(expected)
            assert contains(fit['directed_original_expected_tau_exponent'],model)
            assert contains(fit['fitted_tau_exponent'],model)
            identity=owner.graph.sub(difference,owner.graph.mul(owner.graph.constant(q),expected))
            reported=bg.box.pulse.radius.FunctionRef(owner.graph,fit['exact_scale_exponent_identity'])
            assert not reader.polynomial(reported) and not reader.polynomial(identity)
        assert fit['coefficient_ratio_not_replaced_by_one'] and fit['absolute_width_flags_not_promoted']
        lo,hi=bg.ends(magnitude)
        if fit['actual_magnitude_increase_certified']:assert lo>0
        if fit['actual_magnitude_decrease_certified']:assert hi<0
    for left,right in zip(names,names[1:]):
        comparison=view['time_comparisons'][left+'__to__'+right]
        q=entry[2][right][1].finite_log_time_offset-entry[2][left][1].finite_log_time_offset
        assert q<0 and Fraction(comparison['delta_log_tau'])==q
        tau=owner.graph.sub(entry[2][right][1].forward_coordinates['log_tau'],entry[2][left][1].forward_coordinates['log_tau'])
        assert not reader.polynomial(owner.graph.sub(tau,owner.graph.constant(q)))
        for label,fit in comparison['observable_fits'].items():
            check_fit(fit,entry[3][left][label],entry[3][right][label],q,coefficients[left][label],coefficients[right][label])
            fits_checked+=1
        for label,(a,b) in current.relative.PAIRS.items():
            a='Er_total' if a=='Er' else a
            def quotient(name):
                numerator,denominator=entry[3][name][a],entry[3][name][b]
                return dict(scale=owner.graph.sub(numerator['scale'],denominator['scale']),
                    expected_tau_exponent=None if numerator['expected_tau_exponent'] is None else
                        owner.graph.sub(numerator['expected_tau_exponent'],denominator['expected_tau_exponent']))
            check_fit(comparison['finite_remainder_relative_fits'][label],quotient(left),quotient(right),q,
                coefficients[left][a]/coefficients[left][b],coefficients[right][a]/coefficients[right][b])
            relative_fits_checked+=1
        observations[left+'__to__'+right]=comparison
    assert len(entry[1])>=3 and not view['coefficient_recursion_implemented']
    assert not view['independently_measured_vortex_core_width'] and not view['material_trajectory_or_winding']
    assert not view['all_orders_time_flatness_certified'] and not view['global_or_regional_cone_certified']
    assert view['actual_time_growth_observed']==any(comparison['observable_fits'][label]['actual_magnitude_increase_certified']
        for comparison in entry[4].values() for label in ('Ur','Utheta','Uz','omega_z'))
    assert all(not view[key] for key in current.OPEN)
    invalid=dict(copied_observation=rejected(lambda:owner.report(copy.copy(value))),
        reversed_time_samples=rejected(lambda:owner.evaluate(dict(reversed(list(entry[1].items()))))),
        too_few_samples=rejected(lambda:owner.evaluate(dict(list(entry[1].items())[:2]))))
    first=names[0];sample=entry[1][first];other=entry[1][names[1]]
    mixed=dict(entry[1]);mixed[first]=current.OriginalTimeSample(sample.field,other.axial,sample.background,sample.cone,sample.remainder)
    invalid['mixed_actual_field_packets']=rejected(lambda:owner.evaluate(mixed))
    node=next(iter(entry[5]));saved=copy.deepcopy(owner.graph.nodes[node])
    try:
        owner.graph.nodes[node]['changed_time_scale']=True
        invalid['changed_time_scale_DAG']=rejected(lambda:owner.report(value))
    finally:owner.graph.nodes[node]=saved
    fit=next(iter(entry[4].values()))['observable_fits']['omega_z'];saved=fit['directed_log_magnitude_ratio']
    try:
        fit['directed_log_magnitude_ratio']=owner.ctx.mpf(0)
        invalid['changed_measured_time_fit']=rejected(lambda:owner.report(value))
    finally:fit['directed_log_magnitude_ratio']=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    owner.report(value)
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,actual_multitime_observations=observations,
        actual_sample_count=len(entry[1]),independently_checked_source_products=products_checked,
        independently_checked_point_observable_bounds=bounds_checked,independently_checked_point_cone_bounds=cone_bounds_checked,
        independently_checked_time_fits=fits_checked,independently_checked_remainder_relative_time_fits=relative_fits_checked,
        independent_800_digit_source_product_and_time_fit_bounds=True,
        full_nonzero_Z_lambda_identity_inherited_and_checked=True,
        finite_similarity_path_not_core_width_material_winding_or_coefficient_recursion=True,
        invalid_inputs_rejected=invalid,input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(bg.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_TIME_OBSERVATIONS actual three-time source path and measured fits',flush=True)
    return receipt
