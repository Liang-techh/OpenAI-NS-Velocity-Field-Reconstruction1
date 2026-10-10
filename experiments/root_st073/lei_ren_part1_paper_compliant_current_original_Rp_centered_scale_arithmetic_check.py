"""Independent exact origin arithmetic, finite source boxes and scale widths."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_centered_scale_arithmetic as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def included(outer,inner):
    a,b=current.ends(outer);u,v=current.ends(inner)
    return a<=u<=v<=b


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError,KeyError):return True
    raise AssertionError('Invalid centered source input admitted')


def exact_encoding_identity(singleton,coefficient,record):
    q=Fraction(coefficient)
    sign,mantissa,exponent,bits=singleton[0]
    if not q:
        return record==dict(numerator=0,denominator=1,binary_exponent=0)
    difference=record['binary_exponent']-exponent
    # Only normalization shifts are allocated, never the astronomical origin.
    assert abs(difference)<=abs(q.numerator).bit_length()+q.denominator.bit_length()+mantissa.bit_length()
    actual=Fraction(record['numerator'],record['denominator'])*Fraction(2)**difference
    return actual==(-1 if sign else 1)*mantissa*q


def symbolic(polynomial):
    return s.Add(*(s.Rational(q.numerator,q.denominator)*s.Mul(*(s.Symbol('n'+str(n)) for n in atoms))
                   for atoms,q in polynomial.items()))


@source_precision
def run(owner,deliveries):
    began=time.monotonic()
    raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[k] for k in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name]==digest for name,digest in owner.before.hashes.items())
    independent=MPIntervalContext();independent.dps=500
    c=owner.ctx
    # Direct independent source definitions, not copied new parameter boxes.
    exp40=independent.exp(40);logP=exp40+11
    mu=independent.exp(-4*logP)/1000
    Tw=-60*independent.ln(mu)
    logU0=-independent.mpf(57)/10-exp40/2-(independent.mpf(1)/2+mu)*Tw-mu/2
    finite=dict(mu=mu,logP=logP,Tw=Tw,logU0=logU0,U0=independent.exp(logU0))
    for name,value in finite.items():
        assert included(owner.inclusion[name]['fresh_box'],value),name
        assert included(owner.inclusion[name]['original_box'],owner.inclusion[name]['fresh_box']),name
    # Exact rational origins with very large exponents have bounded mantissas.
    encodings=[]
    for singleton in (owner.exact_singleton,((0,7,10**100,3),)*2,((1,7,-10**100,3),)*2):
        for q in (Fraction(0),Fraction(10),Fraction(-3,8),Fraction(2,3)):
            record=current.exact_scaled_binary_rational(singleton,q)
            assert exact_encoding_identity(singleton,q,record)
            encodings.append(dict(coefficient=str(q),exact_identity=True,
                bounded_mantissa_bits=abs(record['numerator']).bit_length(),
                bounded_exponent_storage_bits=abs(record['binary_exponent']).bit_length()))
    assert rejected(lambda:current.exact_scaled_binary_rational(((0,1,0,1),(0,3,0,2)),Fraction(1)))
    # Finite fixtures make the distinction between arithmetic width and a
    # representation of an exact huge fixed scale observable independently.
    huge=c.mpf(2)**10000
    residual=c.mpf([0,'0.0001'])
    scale=huge+residual
    row=dict(exact_zero_enclosure=False,common_scale_coefficient_enclosure=c.mpf([2,'2.0001']),
        directed_reference_log_scale=scale,exact_reference_log_scale_function=-1,
        ordinary_numeric_materialized=False,requested_relative_width_satisfied=False)
    center=dict(exact_fixed_origin=dict(kind='fixture_exact_binary',value='2^10000'),
        directed_residual_log_scale=residual,residual_log_width_bound=c.mpf(current.ends(residual)[1]),
        exact_centered_scale_identity_proof=-1)
    fixture=current.centered_relative_budget(c,row,Fraction(1,1000),center)
    assert fixture['ordinary_numeric_relative_width_satisfied']
    assert not fixture['ordinary_numeric_delivery_target_satisfied'] and not fixture['ordinary_numeric_materialized']
    assert current.ends(fixture['direct_rounded_reference_log_width_bound'])[0]>1
    physical_ratio=independent.mpf('2.0001')/2*independent.exp(independent.mpf('.0001'))
    # This is an enclosure of an upper width bound, not an enclosure of
    # the exact physical diameter: directed coefficient endpoints may
    # make even its lower endpoint slightly larger than that diameter.
    assert current.ends(physical_ratio-1)[1]<=current.ends(fixture['ordinary_relative_width_bound'])[1]
    broad=dict(center,directed_residual_log_scale=c.mpf([-1,1]),residual_log_width_bound=c.mpf(2))
    assert not current.centered_relative_budget(c,row,Fraction(1,1000),broad)['ordinary_numeric_relative_width_satisfied']
    crossing=dict(row,common_scale_coefficient_enclosure=c.mpf([-1,1]))
    assert current.centered_relative_budget(c,crossing,Fraction(1,1000),center).get('sign') is None
    zeros=dict(row,exact_zero_enclosure=True,ordinary_numeric_materialized=True)
    assert current.centered_relative_budget(c,zeros,Fraction(1,1000),broad)['ordinary_numeric_delivery_target_satisfied']
    observations={};totals=dict(rows=0,exact_zero_rows=0,nonzero_factored_target_rows=0,
        centered_nonzero_rows=0,remaining_nonconstant_origin_rows=0,ordinary_width_target_rows=0,
        ordinary_numeric_target_rows=0,sign_unresolved_rows=0)
    identities=0
    for name,delivery in deliveries.items():
        observed=owner.evaluate(delivery,'1/1000')
        old=(owner.before.evaluate(delivery,'1/1000') if observed['source_coefficient_provider']=='accepted_refined_current_pulse'
             else owner.before.before.evaluate(delivery,'1/1000'))
        reader=owner.reader(delivery)
        for section,components in observed['physical_value_rows'].items():
            for component,rows in components.items():
                for label,row in rows.items():
                    prior=old['physical_value_rows'][section][component][label]
                    assert row['exact_zero_enclosure']==prior['exact_zero_enclosure']
                    if row['exact_zero_enclosure']:continue
                    centered=row['centered_original_log_scale'];budget=row['physical_accuracy']
                    assert row['exact_reference_log_scale_function']==prior['exact_reference_log_scale_function']
                    assert included(prior['directed_reference_log_scale'],row['directed_reference_log_scale'])
                    poly=reader.polynomial(centered['original_full_log_scale_function'])
                    q=poly.get((owner.pivot.node,),Fraction(0));encoding=centered['exact_fixed_origin']
                    record={key:encoding[key] for key in ('numerator','denominator','binary_exponent')}
                    assert exact_encoding_identity(owner.exact_singleton,q,record)
                    source=reader.polynomial(encoding['source_function'])
                    residual=reader.polynomial(centered['exact_residual_log_scale_function'])
                    assert s.expand(symbolic(source)+symbolic(residual)-symbolic(poly))==0
                    remaining=any(owner.pivot.node in atoms for atoms in residual)
                    assert remaining==centered['residual_still_contains_nonconstant_logC_products']
                    # Independent polynomial interval arithmetic over the same
                    # actual original atomic source boxes, including products
                    # still involving logC and a nonconstant coefficient.
                    independent_residual=independent.mpf(0)
                    for atoms,coefficient in residual.items():
                        term=independent.mpf(coefficient.numerator)/coefficient.denominator
                        for atom in atoms:term*=independent.mpf(current.ends(reader.at(atom)))
                        independent_residual+=term
                    assert included(centered['directed_residual_log_scale'],independent_residual)
                    a,b=current.ends(centered['directed_residual_log_scale'])
                    exact_width=independent.mpf(b)-independent.mpf(a)
                    assert included(centered['residual_log_width_bound'],exact_width)
                    assert included(row['common_scale_coefficient_enclosure'],
                        independent.mpf(current.ends(row['common_scale_coefficient_enclosure'])))
                    if budget.get('sign'):
                        a,b=current.ends(row['common_scale_coefficient_enclosure'])
                        magnitude_a=abs(independent.mpf(a));magnitude_b=abs(independent.mpf(b))
                        lo=min(current.ends(magnitude_a)[0],current.ends(magnitude_b)[0])
                        hi=max(current.ends(magnitude_a)[1],current.ends(magnitude_b)[1])
                        independent_ratio=exact_width+independent.ln(independent.mpf(hi)/independent.mpf(lo))
                        assert included(budget['log_ordinary_magnitude_ratio_bound'],independent_ratio)
                        assert not budget['ordinary_numeric_delivery_target_satisfied']
                    assert encoding['zero_error_variation'] and encoding['no_exponent_sized_integer_shift']
                    assert budget['exact_fixed_origin_is_not_a_numeric_scale_approximation']
                    assert budget['normalized_error_does_not_replace_absolute_physical_error']
                    identities+=1
        for key,value in observed['delivery_counts'].items():totals[key]+=value
        packet=owner.velocity_pressure(delivery,'1/1000')
        assert set(packet['values'])=={'u','v','w','p'}
        assert not packet['full_certified_physical_accuracy'] and not packet['unrestricted_physical_point_API']
        assert not any(packet[key] for key in current.OPEN)
        observations[name]=dict(counts=observed['delivery_counts'],base_accuracy={component:{
            key:value for key,value in row['physical_accuracy'].items()
            if key not in ('actual_common_scale_source_ledger','signed_coefficient_and_ratio_variation_ledger')}
            for component,row in packet['values'].items()},
            old_base_reference_widths={component:prior['physical_accuracy'].get('reference_log_width_bound')
                for component,key in (('u','ux'),('v','uy'),('w','uz'),('p','p'))
                for prior in [old['physical_value_rows']['Cartesian_spatial_rows'][key]['x0_y0_z0']]},
            history_partition_scope=observed['coefficient_history_partition'])
    first=next(iter(deliveries.values()));reader=owner.reader(first)
    log_function=owner.before.velocity_pressure(first,'1/1000')['values']['u']['exact_reference_log_scale_function']
    original_center=owner.center(reader,log_function)
    # Memo injection is ignored; only sealed source bindings/aliases matter.
    reader.memo[original_center['exact_residual_log_scale_function']]=c.mpf(0)
    repeated=owner.center(reader,log_function)
    assert repeated['directed_residual_log_scale']._mpi_==original_center['directed_residual_log_scale']._mpi_
    invalid=dict(copied_delivery=rejected(lambda:owner.evaluate(copy.copy(first))),
        invalid_target=rejected(lambda:owner.evaluate(first,'0')),
        copied_reader=rejected(lambda:owner.center(copy.copy(reader),log_function)),
        foreign_graph_function=rejected(lambda:owner.center(reader,
            current.refined.box.pulse.radius.FunctionRef(object(),log_function))))
    saved=reader.bindings[owner.amplitude.mu.node]
    try:
        reader.bindings[owner.amplitude.mu.node]=c.mpf(1)
        invalid['reader_parameter_oracle']=rejected(lambda:owner.center(reader,log_function))
    finally:reader.bindings[owner.amplitude.mu.node]=saved
    saved=dict(reader.aliases)
    try:
        reader.aliases[log_function]=owner.graph.zero.node
        invalid['reader_alias_oracle']=rejected(lambda:owner.center(reader,log_function))
    finally:reader.aliases.clear();reader.aliases.update(saved)
    saved=owner.finite_bindings[owner.amplitude.mu.node]
    try:
        owner.finite_bindings[owner.amplitude.mu.node]=c.mpf(1)
        invalid['owner_parameter_substitution']=rejected(owner.assert_graph)
    finally:owner.finite_bindings[owner.amplitude.mu.node]=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME)
    hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,actual_row_counts=totals,
        actual_exact_centered_identities=identities,actual_centered_observations=observations,
        independent_500_digit_elementary_parameter_inclusion=True,
        independent_exact_rational_binary_origins=encodings,
        finite_physical_ratio_fixture_included=True,
        exact_origin_and_residual_width_not_numeric_scale_materialization=True,
        nonconstant_logC_products_retained_in_residual=True,
        invalid_inputs_rejected=invalid,caller_memo_oracle_ignored=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_CENTERED_SCALE_ARITHMETIC exact fixed origins and honest residual error',flush=True)
    return receipt
