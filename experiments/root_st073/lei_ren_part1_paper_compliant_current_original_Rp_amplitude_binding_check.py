"""Independent reflection, source parameter and current U amplitude checks."""
import copy
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_amplitude_binding as current
import lei_ren_part1_paper_compliant_current_original_Rp_native_constants_check as native_check
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError):return True
    raise AssertionError('Invalid current amplitude input was admitted')


@source_precision
def run(owner=None,deliveries=None):
    began=time.monotonic();raw=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    assert not any(raw[key] for key in current.GATES+current.OPEN)
    owner=owner if owner is not None else current.CurrentOriginalRpAmplitudeBinding(require_checked=False)
    assert not owner.acceptance_loaded and all(owner.assert_graph().values())
    observed=current.correlated.report(owner.report())
    for key,value in observed.items():assert raw[key]==value,key
    q,constants=owner.exact.current_constants();_,kernels=owner.exact.original_native_functions(q)
    # The actual source AST supplies its factors; this checker does not read
    # the new compact expression as the definition of the original U.
    projector,_,coefficients,stages=native_check.source_frontend(owner.exact,q,kernels)
    assert current.zero(coefficients['U']-constants['U'])
    initial=s.exp(s.Rational(1,10)-s.Rational(3,5)*kernels['J1'])
    axial=s.exp(-(s.exp(40)+10)/2)
    transition=s.exp(-s.Rational(1,2)-q.mu*kernels['J_transition'])
    Tw=q.at(owner.exact.frame.functions['Tw'])
    power=s.exp(-(s.Rational(1,2)+q.mu)*Tw)
    assert current.zero(coefficients['U']-initial*axial*transition*power)
    a,b=s.symbols('a b',positive=True)
    assert s.cancel(a/(a+b)+b/(a+b)-1)==0
    endpoint_subs={kernels[key]:s.Rational(1,2) for key in ('J1','J_transition')}
    independent_log=s.Rational(1,10)-s.Rational(3,5)*s.Rational(1,2)
    independent_log-=(s.exp(40)+10)/2+s.Rational(1,2)+q.mu/2+(s.Rational(1,2)+q.mu)*Tw
    assert current.zero(coefficients['U'].subs(endpoint_subs)-s.exp(independent_log))
    proof=current.sigma_endpoint_proofs(owner.exact,q,kernels)
    assert proof==owner.sigma_proofs
    predicates=current.exact.frame.bridge.leading_check.source_predicates()
    # Read exact source parameter nodes without the usual formal mu/logP
    # bindings; equality is to their defining expressions, never box ends.
    symbolic=current.exact.frame.bridge.leading_check.Interpreter(owner.exact.frame.bridge.leading)
    symbolic.g=owner.graph
    symbolic.bindings.pop(owner.radius.parameters['mu'].node)
    symbolic.bindings.pop(owner.radius.parameters['logP'].node)
    exact_logP=s.exp(40)+11;exact_mu=s.exp(-4*exact_logP)/1000
    assert current.zero(symbolic.at(owner.logP)-exact_logP)
    assert current.zero(symbolic.at(owner.mu)-exact_mu)
    assert s.simplify(symbolic.at(owner.Tw)+60*s.log(exact_mu))==0
    # Independent directed evaluation at 400 digits, using the separated
    # source decay factors rather than the new graph compact-log route.
    c=MPIntervalContext();c.dps=400
    ep=c.exp(c.mpf(40));lp=ep+11;lm=-c.ln(c.mpf(1000))-4*lp
    mu=c.exp(lm);tw=-60*lm
    source_log=c.mpf(1)/10-(c.mpf(3)/5)/2-(ep+10)/2-c.mpf(1)/2-mu/2-(c.mpf(1)/2+mu)*tw
    source_U=c.exp(c.mpf(1)/10-(c.mpf(3)/5)/2)*c.exp(-(ep+10)/2)
    source_U*=c.exp(-c.mpf(1)/2-mu/2)*c.exp(-(c.mpf(1)/2+mu)*tw)
    assert current.contains(owner.logU0_box,source_log)
    assert current.contains(owner.U0_box,source_U)
    assert current.contains(owner.inlet.constants['U'],source_U)
    assert current.contains(owner.selected.constants['U'],source_U)
    assert current.contains(owner.logP_box,lp) and current.contains(owner.Tw_box,tw)
    reader=owner.reader()
    assert current.contains(reader.at(owner.raw.U0),source_U)
    assert current.contains(reader.at(owner.raw.logU0),source_log)
    assert current.contains(reader.at(owner.logP),lp)
    outputs={};row_count=group_count=0
    if not deliveries:raise ValueError('Actual live correlated source delivery is required for the scale-reader gate')
    for name,delivery in deliveries.items():
        out=owner.physical_log_scale_rows(delivery);outputs[name]=out
        actual=owner.before.physical_rows(delivery);numeric=owner.reader(delivery)
        for section,components in out['evaluated_physical_log_scale_rows'].items():
            for component,rows in components.items():
                originals=actual[section][component]
                originals=originals if section=='Cartesian_spatial_rows' else {'dt':originals}
                assert set(rows)==set(originals)
                for label,row in rows.items():
                    source=originals[label];assert row['physical_derivative']==source.derivative
                    groups=source.groups();assert len(row['groups'])==len(groups);row_count+=1
                    for report,group in zip(row['groups'],groups):
                        exact=owner.graph.add(*[ref for _,ref in group['log_scale_parts']])
                        assert report['exact_positive_scale_log_function']==exact.node
                        assert report['signed_source_coefficient'].coefficients==group['signed_coefficient'].coefficients
                        assert report['source_scale_units']==group['scale_units']
                        assert report['source_scale_powers']==[dict(numerator=v.numerator,denominator=v.denominator)
                            for v in group['exact_source_scale_powers']]
                        assert report['directed_positive_scale_log']._mpi_==numeric.at(exact)._mpi_
                        assert all(mp.isfinite(v) for v in current.ends(report['directed_positive_scale_log']))
                        group_count+=1
        assert out['absolute_physical_values_and_accuracy_not_materialized']
        assert not any(out[key] for key in current.GATES+current.OPEN)
    bad={'copied_delivery':rejected(lambda:owner.reader(copy.copy(next(iter(deliveries.values()))))),
         'caller_numeric_bound':rejected(lambda:owner.reader({'logU0':[0,0]})),
         'unproved_absolute_radius':rejected(lambda:reader.at(owner.graph.unary('exp',owner.radius.logRp)))}
    saved=owner.logU0
    try:
        owner.logU0=owner.graph.zero
        bad['mutated_compact_alias']=rejected(lambda:owner.reader())
    finally:owner.logU0=saved
    saved_box=owner.U0_box
    try:
        owner.U0_box=owner.amplitude_ctx.mpf(1)
        bad['mutated_U0_bound']=rejected(lambda:owner.reader())
    finally:owner.U0_box=saved_box
    assert all(owner.assert_graph().values()) and all(bad.values())
    result=dict(all_passed=True,source_family=owner.family_record,
        actual_source_AST_amplitude_factor_identity=True,actual_source_stages=stages,
        reflected_flat_sigma_endpoint_integral_proofs=proof,
        independent_source_defining_predicates=predicates,
        actual_exact_logP_mu_Tw_parameter_definitions_checked=True,
        independent_400_digit_separated_source_U_and_logU_included=True,
        full_existing_inlet_and_selected_U_boxes_include_actual_current_U0=True,
        actual_current_physical_operator_rows_with_directed_log_scales=row_count,
        actual_signed_physical_scale_groups=group_count,
        actual_live_source_scale_rows=outputs,invalid_inputs_rejected=bad,
        original_source_owners_and_full_boxes_unmutated=True,
        P0_inclusion_and_distinct_scale_sum_physical_accuracy_remain_open=True,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,**projector.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(result),indent=2)+'\n',
        encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_U0_SOURCE_BINDING source amplitude and physical log scales',flush=True)
    return result


if __name__=='__main__':run()
