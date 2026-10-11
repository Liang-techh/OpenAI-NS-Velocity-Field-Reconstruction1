"""Independent source hierarchy, mixed-product bounds and physical widths."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_source_product_arithmetic as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError,KeyError):return True
    raise AssertionError('Invalid source-product input admitted')


def expression(graph,node):
    node=node.node if hasattr(node,'node') else node
    n=graph.nodes[node];op=n['operation']
    if op=='exact_rational':return s.Rational(n['numerator'],n['denominator'])
    if op=='sum':return s.Add(*(expression(graph,v) for v in n['arguments']))
    if op=='product':return s.Mul(*(expression(graph,v) for v in n['arguments']))
    if op=='negative':return -expression(graph,n['argument'])
    if op=='positive_quotient':return expression(graph,n['numerator'])/expression(graph,n['denominator'])
    if op=='analytic_unary' and n['name'] in ('log','exp'):
        return getattr(s,n['name'])(expression(graph,n['argument']))
    raise ValueError('Closed original elementary source expression required')


def symbolic(poly):
    return s.Add(*(s.Rational(q.numerator,q.denominator)*s.Mul(*(s.Symbol('n'+str(v)) for v in atoms))
        for atoms,q in poly.items()))


@source_precision
def run(owner,deliveries):
    began=time.monotonic()
    raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name]==digest for name,digest in owner.before.hashes.items())
    c=MPIntervalContext();c.dps=800
    ep=c.exp(40);P=ep+11;logmu=c.ln(c.mpf(1)/1000)-4*P
    mu=c.exp(-4*P)/1000;ld=-4*P-30;delta=c.exp(ld);Tw=-60*logmu
    logU=-c.mpf(57)/10-ep/2-(c.mpf(1)/2+mu)*Tw-mu/2
    independent=dict(mu=mu,logP=P,Tw=Tw,logU0=logU,delta=delta,logmu=logmu,
        logdelta=ld,U0=c.exp(logU))
    expectedP=s.exp(40)+11
    assert s.simplify(expression(owner.graph,owner.amplitude.logP)-expectedP)==0
    assert s.simplify(expression(owner.graph,owner.delta)-s.exp(-4*expectedP-30))==0
    assert s.simplify(expression(owner.graph,owner.amplitude.mu)-s.exp(-4*expectedP)/1000)==0
    a=s.Symbol('a',real=True)
    assert s.simplify(s.exp(-4*a-30)/(s.exp(-4*a)/1000)-1000*s.exp(-30))==0
    assert current.ends(ld)[1]<current.ends(-200*c.ln(10))[0]
    assert owner.native_definition['delta']=='min(1e-200,exp(-4logPstar-30))'
    for name,value in independent.items():
        assert current.contains(owner.inclusion[name]['fresh_box'],value),name
        assert current.contains(owner.inclusion[name]['original_box'],owner.inclusion[name]['fresh_box']),name
    ratio=c.mpf(1000)*c.exp(-30)
    actual_ratio=owner.ctx.mpf(owner.finite_bindings[owner.delta.node])/owner.finite_bindings[owner.amplitude.mu.node]
    assert current.contains(actual_ratio,ratio)
    totals={};observations={};identities=0;recovered=0
    for name,delivery in deliveries.items():
        observed=owner.evaluate(delivery,'1/1000');old=owner.before.evaluate(delivery,'1/1000')
        reader=owner.reader(delivery);original=owner.before.before.reader(delivery)
        for node,value in owner.finite_bindings.items():
            assert current.contains(original.at(node),value),node
        bindings=dict(original.bindings)
        for key,ref in owner.refs.items():bindings[ref.node]=independent[key]
        bindings[owner.amplitude.raw.U0.node]=independent['U0']
        independent_reader=current.correlated.CorrelatedGraphBounds(
            current.correlated.locator.CancelledGraphBounds(owner.graph,c,bindings,original.allowed_exponentials),original.aliases)
        product=owner.ctx.mpf(reader.bindings[owner.pivot.node])*reader.at(owner.delta)
        independent_product=c.mpf(current.ends(reader.bindings[owner.pivot.node]))*delta
        assert current.contains(product,independent_product)
        assert current.ends(product)[0]>0 and current.ends(product)[0]!=current.ends(product)[1]
        base_flags={}
        for section,components in observed['physical_value_rows'].items():
            for component,rows in components.items():
                for label,row in rows.items():
                    prior=old['physical_value_rows'][section][component][label]
                    assert row['exact_reference_log_scale_function']==prior['exact_reference_log_scale_function']
                    assert current.contains(prior['directed_reference_log_scale'],row['directed_reference_log_scale'])
                    assert current.contains(prior['common_scale_coefficient_enclosure'],row['common_scale_coefficient_enclosure'])
                    assert row['physical_derivative']==prior['physical_derivative']
                    new_center=row['centered_original_log_scale'];old_center=prior['centered_original_log_scale']
                    assert new_center['exact_fixed_origin']==old_center['exact_fixed_origin']
                    poly=reader.polynomial(new_center['original_full_log_scale_function'])
                    fixed=reader.polynomial(new_center['exact_fixed_origin']['source_function'])
                    residual=reader.polynomial(new_center['exact_residual_log_scale_function'])
                    assert s.expand(symbolic(poly)-symbolic(fixed)-symbolic(residual))==0
                    independent_residual=c.mpf(0)
                    for atoms,q in residual.items():
                        term=c.mpf(q.numerator)/q.denominator
                        for atom in atoms:term*=independent_reader.at(atom)
                        independent_residual+=term
                    assert current.contains(new_center['directed_residual_log_scale'],independent_residual)
                    assert current.contains(old_center['directed_residual_log_scale'],new_center['directed_residual_log_scale'])
                    l,h=current.ends(new_center['directed_residual_log_scale']);width=c.mpf(h)-c.mpf(l)
                    assert current.contains(new_center['residual_log_width_bound'],width)
                    lo,hi=current.ends(row['common_scale_coefficient_enclosure'])
                    # Keep endpoint magnitudes exact before independent arithmetic.
                    magnitudes=[mp.make_mpf((0,*v._mpf_[1:])) for v in (lo,hi)]
                    coeff_ratio=c.mpf(max(magnitudes))/c.mpf(min(magnitudes))
                    physical_log_ratio=width+c.ln(coeff_ratio)
                    budget=row['physical_accuracy'];allowed=c.ln(1+c.mpf(1)/1000)
                    assert current.contains(budget['log_ordinary_magnitude_ratio_bound'],physical_log_ratio)
                    if budget['ordinary_numeric_relative_width_satisfied']:
                        assert current.ends(physical_log_ratio)[1]<=current.ends(allowed)[0]
                        assert current.ends(coeff_ratio*c.exp(width)-1)[1]<=current.ends(budget['ordinary_relative_width_bound'])[1]
                        recovered+=int(not prior['physical_accuracy']['ordinary_numeric_relative_width_satisfied'])
                    assert budget['sign']==prior['physical_accuracy']['sign']
                    assert not budget['ordinary_numeric_delivery_target_satisfied'] and not row['ordinary_numeric_materialized']
                    assert new_center['residual_still_contains_nonconstant_logC_products']==old_center['residual_still_contains_nonconstant_logC_products']
                    identities+=1
                    if section=='Cartesian_spatial_rows' and label=='x0_y0_z0':
                        base_flags[component]=dict(ordinary_relative_width_satisfied=budget['ordinary_numeric_relative_width_satisfied'],
                            coefficient_relative_width_upper=budget['coefficient_relative_width_upper'],
                            ordinary_relative_width_bound=budget.get('ordinary_relative_width_bound'),
                            centered_residual_width=new_center['residual_log_width_bound'],
                            actual_numeric_materialization=False)
        assert base_flags['uz']['ordinary_relative_width_satisfied'] and base_flags['p']['ordinary_relative_width_satisfied']
        assert not base_flags['ux']['ordinary_relative_width_satisfied'] and not base_flags['uy']['ordinary_relative_width_satisfied']
        for key,value in observed['delivery_counts'].items():totals[key]=totals.get(key,0)+value
        observations[name]=dict(counts=observed['delivery_counts'],old_counts=old['delivery_counts'],
            base_width_flags=base_flags,delta_logC_product_enclosure=product,
            mixed_product_uncertainty_positive=True,same_original_physical_scale_functions=True)
    first=next(iter(deliveries.values()))
    invalid=dict(copied_delivery=rejected(lambda:owner.evaluate(copy.copy(first))),
        invalid_target=rejected(lambda:owner.evaluate(first,'0')))
    reader=owner.reader(first)
    # Existing centered API must ignore caller memo/polynomial oracles.
    scale=owner.before.velocity_pressure(first,'1/1000')['values']['w']['exact_reference_log_scale_function']
    clean=owner.center(reader,scale)
    reader.memo[clean['exact_residual_log_scale_function']]=owner.ctx.mpf(0)
    reader.polys[scale]={}
    assert clean['directed_residual_log_scale']._mpi_==owner.center(reader,scale)['directed_residual_log_scale']._mpi_
    reader.bindings[owner.delta.node]=owner.ctx.mpf(0)
    invalid['reader_parameter_substitution']=rejected(lambda:owner.center(reader,scale))
    saved=owner.finite_bindings.pop(owner.delta.node)
    try:invalid['missing_delta_binding']=rejected(owner.assert_graph)
    finally:owner.finite_bindings[owner.delta.node]=saved
    saved=copy.deepcopy(owner.graph.nodes[owner.delta.node])
    try:
        owner.graph.nodes[owner.delta.node]['argument']=owner.amplitude.mu.node
        invalid['delta_child_DAG_substitution']=rejected(owner.assert_graph)
    finally:owner.graph.nodes[owner.delta.node].clear();owner.graph.nodes[owner.delta.node].update(saved)
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME)
    hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,actual_row_counts=totals,
        actual_source_product_observations=observations,actual_scale_identities_checked=identities,
        recovered_ordinary_relative_width_rows=recovered,
        independent_800_digit_original_source_parameter_definitions=True,
        independent_exact_mu_delta_ratio_and_strict_min_branch=True,
        independent_residual_polynomials_and_physical_error_budgets=True,
        original_fixed_origins_scales_coordinates_pressure_and_signed_sources_retained=True,
        mixed_source_products_keep_nonzero_directed_uncertainty=True,
        no_physical_exponential_guard_or_materialization_promoted=True,
        invalid_inputs_rejected=invalid,input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_SOURCE_PRODUCT_ARITHMETIC original delta products and relative scale widths',flush=True)
    return receipt
