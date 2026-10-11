"""Independent source products, complete physical scales and relative bounds."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_remainder_ratios as current
import lei_ren_part1_paper_compliant_current_original_Rp_background_stress_check as stress_check
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

bg=current.bg


def original_reader(owner,delivery,c):
    issued=owner.product.reader(delivery)
    return bg.correlated.CorrelatedGraphBounds(bg.correlated.locator.CancelledGraphBounds(
        owner.graph,c,dict(issued.bindings),issued.allowed_exponentials),dict(issued.aliases))


def independent_row(owner,delivery,row,c,*,source=None,request=None,reader=None):
    # A checker may reuse one already validated live source/reader within a
    # delivery. It validates the packet again after the independent arithmetic.
    if source is None:source=owner.product.before.source(delivery)
    if request is None:_,_,request=owner.product.amplitude.before._validate_delivery(delivery)
    if reader is None:reader=original_reader(owner,delivery,c)
    angle=c.mpf(request.theta.numerator)/request.theta.denominator
    bindings={bg.Z:c.mpf(request.Z.numerator)/request.Z.denominator,
        bg.DELTA:c.mpf(reader.bindings[owner.product.delta.node]),
        bg.MU:c.mpf(reader.bindings[owner.product.amplitude.mu.node]),bg.CS:c.cos(angle),bg.SN:c.sin(angle)}
    groups={};products=0
    for item in row['actual_source_product_ledger']:
        coefficient=stress_check.numeric(s.sympify(item['exact_operator_expression']),c,bindings)
        units={'R':Fraction(item['radial_power']['numerator'],item['radial_power']['denominator'])}
        product=c.mpf(1)
        for atom in item['actual_shared_source_product']:
            label=atom['source'];k,n=atom['ordinary_derivative'];power=atom['power']
            primitive=source['log_radius_mixed_rows'][label]['y%d_Z%d'%(k,n)]
            assert atom['actual_source_units']==primitive.source_units
            assert atom['actual_source_powers']==[bg.signed.rational_record(q) for q in primitive.powers]
            product*=c.mpf(primitive.coefficients[0])**power
            for unit,q in zip(primitive.source_units,primitive.powers):units[unit]=units.get(unit,Fraction(0))+power*q
        assert bg.axial.differential.product.contains(item['signed_product_coefficient'],coefficient*product)
        log=bg.box.pulse.radius.FunctionRef(owner.graph,item['exact_product_log_scale'])
        key=(tuple(sorted(reader.polynomial(log).items())),tuple(sorted(units.items())),item['lambda_exponent'])
        group=groups.setdefault(key,dict(log=log,coefficient=c.mpf(0)))
        group['coefficient']+=coefficient*product;products+=1
    if row['exact_zero_enclosure']:
        assert not groups and bg.ends(row['ordinary_numeric_enclosure'])==(0,0)
        return c.mpf(0),products
    reference=bg.box.pulse.radius.FunctionRef(owner.graph,row['exact_reference_log_scale_function'])
    total=c.mpf(0);failed=False
    for group in groups.values():
        difference=owner.graph.sub(group['log'],reference)
        if not reader.polynomial(difference):ratio=c.mpf(1)
        else:
            lo,hi=bg.ends(reader.at(difference))
            if hi>1000:failed=True;continue
            if hi<=-1000:ratio=c.mpf([0,bg.ends(c.exp(-1000))[1]])
            elif lo<-1000:ratio=c.mpf([0,bg.ends(c.exp(c.mpf(hi)))[1]])
            else:ratio=c.exp(reader.at(difference))
        total+=group['coefficient']*ratio
    if failed:
        assert row['common_scale_coefficient_enclosure'] is None
        return None,products
    assert bg.axial.differential.product.contains(row['common_scale_coefficient_enclosure'],total)
    return total,products


def independent_ratio(c,coefficient,log_bound):
    lo,hi=bg.ends(log_bound)
    if hi>1000:return None
    if hi<=-1000:scale=c.mpf([0,bg.ends(c.exp(-1000))[1]])
    elif lo<-1000:scale=c.mpf([0,bg.ends(c.exp(c.mpf(hi)))[1]])
    else:scale=c.exp(log_bound)
    return coefficient*scale


@source_precision
def run(owner,values):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not raw[current.GATE] and not any(raw[key] for key in current.OPEN)
    assert raw['source_family']==owner.family_record and raw['original_ratio_pairs']=={k:list(v) for k,v in current.PAIRS.items()}
    assert raw['input_hashes']==owner.hashes
    assert raw['input_hashes'][Path(current.__file__).name]==current.sha(Path(current.__file__).name)
    assert raw['actual_original_finite_remainder_ratios']==bg.correlated.report({name:owner.report(value) for name,value in values.items()})
    assert all(owner.assert_graph().values())
    c=MPIntervalContext();c.dps=800;observations={};products_checked=ratios_checked=0;nonzero_z_seen=False
    for name,value in values.items():
        report=owner.report(value);entry=owner._fields[id(value)]
        field=owner.background._fields[id(entry[1])][1]
        record=owner.background.differential._require(field);delivery=record['delivery']
        _,_,request=owner.product.amplitude.before._validate_delivery(delivery)
        nonzero_z_seen=nonzero_z_seen or request.Z!=0
        reader=original_reader(owner,delivery,c);coordinates=request.forward_coordinates
        identity=owner.graph.add(owner.graph.mul(owner.graph.constant(2),coordinates['loglambda']),
            owner.graph.neg(coordinates['log_tau']),owner.graph.unary('log',owner.graph.sub(owner.graph.one,
                owner.graph.mul(coordinates['Z'],coordinates['Z']))))
        assert not reader.polynomial(identity)
        independent={};source=owner.product.before.source(delivery)
        for label,row in entry[2].items():
            coefficient,count=independent_row(owner,delivery,row,c,source=source,request=request,reader=reader)
            assert coefficient is not None
            independent[label]=coefficient;products_checked+=count
        ratio_observations={}
        for label,(numerator,denominator) in current.PAIRS.items():
            result=report['source_correlated_ratios'][label]
            a=owner.arithmetic._require(entry[3][numerator]);b=owner.arithmetic._require(entry[3][denominator])
            assert a['delivery'] is b['delivery'] is delivery and a['unit']==b['unit']=='acceleration'
            assert b['sign']==1 and bg.ends(independent[denominator])[0]>0
            difference=owner.graph.sub(a['scale'],b['scale'])
            reported=bg.box.pulse.radius.FunctionRef(owner.graph,result['exact_scale_log_difference'])
            assert not reader.polynomial(owner.graph.sub(reported,difference))
            coefficient=independent[numerator]/independent[denominator]
            assert bg.axial.differential.product.contains(result['signed_coefficient_ratio_enclosure'],coefficient)
            bound=reader.at(difference)
            assert bg.axial.differential.product.contains(result['directed_scale_log_difference'],bound)
            numeric=independent_ratio(c,coefficient,bound)
            scale=independent_ratio(c,c.mpf(1),bound)
            if scale is None:assert result['directed_scale_ratio_enclosure'] is None
            else:assert bg.axial.differential.product.contains(result['directed_scale_ratio_enclosure'],scale)
            if numeric is None:assert result['ordinary_numeric_ratio_enclosure'] is None
            else:assert bg.axial.differential.product.contains(result['ordinary_numeric_ratio_enclosure'],numeric)
            assert result['sign']==a['sign']*b['sign'] and result['exact_ratio_is_nonzero']
            assert result['absolute_physical_accuracy_not_promoted']
            ratio_observations[label]=dict(sign=result['sign'],ratio_method=result['ratio_method'],
                signed_coefficient_ratio=result['signed_coefficient_ratio_enclosure'],
                directed_log_scale_difference=result['directed_scale_log_difference'],
                ordinary_numeric_ratio=result['ordinary_numeric_ratio_enclosure'])
            ratios_checked+=1
        assert report['flatness_not_certified'] and not report['all_orders_time_flatness_certified']
        assert report['actual_complete_lambda_at_nonzero_Z_retained']==(request.Z!=0)
        assert all(not report[key] for key in current.OPEN)
        observations[name]=dict(chart=request.chart,native=str(request.native),Z=str(request.Z),theta=str(request.theta),
            finite_log_time_offset=str(request.finite_log_time_offset),ratios=ratio_observations,
            finite_sector_point_only=True,flatness_not_certified=True)
    value=next(iter(values.values()));entry=owner._fields[id(value)]
    assert nonzero_z_seen
    invalid=dict(copied_ratio_packet=stress_check.rejected(lambda:owner.report(copy.copy(value))),
        invalid_target=stress_check.rejected(lambda:owner.evaluate(entry[1],'2')))
    node=next(iter(entry[6]));saved=copy.deepcopy(owner.graph.nodes[node])
    try:
        owner.graph.nodes[node]['changed_finite_remainder_ratio']=True
        invalid['changed_ratio_defining_DAG']=stress_check.rejected(lambda:owner.report(value))
    finally:owner.graph.nodes[node]=saved
    issued=entry[3]['Etheta'];record=owner.arithmetic._values[id(issued)][1];saved=record['coefficient_tuple']
    try:
        record['coefficient_tuple']=owner.ctx.mpf(1)._mpi_
        invalid['changed_issued_finite_remainder_coefficient']=stress_check.rejected(lambda:owner.report(value))
    finally:record['coefficient_tuple']=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    for value in values.values():owner.report(value)
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,original_ratio_pairs={k:list(v) for k,v in current.PAIRS.items()},
        actual_finite_remainder_relative_observations=observations,independently_checked_shared_source_products=products_checked,
        independently_checked_relative_bounds=ratios_checked,independent_800_digit_source_product_and_ratio_bounds=True,
        exact_nonzero_Z_lambda_in_all_physical_scales=True,finite_sector_not_all_orders_time_flat=True,
        invalid_inputs_rejected=invalid,input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(bg.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_REMAINDER_RATIOS finite physical-sector relative bounds',flush=True)
    return receipt
