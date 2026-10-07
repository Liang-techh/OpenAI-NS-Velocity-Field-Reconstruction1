"""Focused current-source attachment, analytic envelope and open-inlet checks."""
import json
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_inner_relaxed_inputs as source


def source_bound_fixtures():
    """Independent exact-rational envelope cases; never current hb values.

    Exact arithmetic includes equality at a=hb/(2K), so a rounded
    boundary product cannot create a false failure of the analytic bound.
    """
    count=0
    for K in (s.Integer(10)**6,s.Integer(10)**9):
        hb=s.Rational(1,10**1000);lower=hb/(2*K)
        for Db in (1/(2*K),s.Rational(7,2),K):
            for Eb in (-K,s.Integer(0),K):
                H=(Db*Db+Eb*Eb)/Db
                if H>K**10:continue
                for sig in (s.Integer(0),s.Rational(37,100),s.Integer(1)):
                    chi=1-(1-hb)*sig
                    phases=((chi*Db,-chi*Eb),(hb*Db,-hb*Eb*(1-sig)),
                        (hb*Db*(1-sig)+s.Rational(4,5)*sig,s.Integer(0)),(s.Rational(4,5),s.Integer(0)))
                    for a,b in phases:
                        if a<lower or a+b*b/a>K**10 or abs(b/a)>4*K*K:
                            raise ArithmeticError('Independent analytic a/kappa/t0 envelope failed')
                        count+=3
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed inner source attachment: '+name)
    owner=source.CurrentInnerRelaxedInputs()
    theorem=source.exact_switch_source_theorem()
    if source.inputs.packets.encode(theorem)!=data['exact_original_micro_power_and_H0_source_theorem']:
        raise ValueError('Actual original switch/power source identities changed')
    if source.inputs.packets.encode(owner.proof)!=data['current_whole_inner_relaxed_input_and_bounds']:
        raise ValueError('Current analytic inner source bounds changed')
    if owner.conditions!=data['current_actual_inner_source_function_attachment_conditions']:
        raise ValueError('Current callable/source identity attachment changed')
    query_count=0
    for chart in source.CHARTS:
        coordinate=('.01','1') if chart=='bridge_first' else source.inputs.packets.DOMAINS[chart]
        query=owner.query(chart,coordinate)
        if not query[source.GATE] or not query['relaxed_input_certified_for_every_source_point_in_query']:
            raise ArithmeticError('Actual analytic inner chart source gate missing')
        if source.inputs.packets.encode(query)!=data['analytic_source_subbox_examples'][chart]:
            raise ValueError('Actual analytic source subbox changed')
        query_count+=1
    for coordinate,key in (((0,1),'closed_first_with_separate_zero'),(0,'exact_zero_Ra')):
        query=owner.query('bridge_first',coordinate)
        if query[source.GATE] or not query['excluded_Ra_stress_free_endpoint_may_be_in_query']:
            raise ArithmeticError('Zero core inlet was promoted to strict relaxed input')
        if source.inputs.packets.encode(query)!=data['analytic_source_subbox_examples'][key]:
            raise ValueError('Separate core zero/source query changed')
        query_count+=1
    bad=(('core',1),('bridge_first',-1),('bridge_second',0),('switch_power',2),('bridge_macro','inf'))
    for chart,coordinate in bad:
        try:owner.query(chart,coordinate)
        except ValueError:pass
        else:raise ArithmeticError('Invalid original inner chart/query admitted')
    try:owner.query('bridge_first','.5',Z=2)
    except ValueError:pass
    else:raise ArithmeticError('Invalid whole axial source domain admitted')
    if any(data[k] or owner.proof[k] for k in source.OPEN):raise ArithmeticError('Scoped inner source gate promoted global targets')
    result=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        source_family=owner.family,current_actual_source_identity_conditions=len(owner.conditions),
        exact_original_micro_power_H0_identities=len(theorem['exact_shear_source_identities']),
        exact_original_control_and_power_AST_bindings=len(theorem['original_switch_and_postpower_AST_bindings']),
        positive_current_analytic_margins=len(owner.proof['positive_directed_margins']),
        independent_a_kappa_t0_source_envelope_inequalities=source_bound_fixtures(),
        source_subboxes_checked=query_count,invalid_original_source_queries_rejected=len(bad)+1,
        source_width_and_a_lower_bound_not_materialized=True,
        zero_core_inlet_case_not_promoted=True,strict_completed_tensor_cone_new_regions_admitted=0,
        source_ancestor_constructors_called=False,p1_p2_whole_inner_norm_bounds_certified=False,
        scope='Whole current original analytic open inner exit through R110; relaxed-only, actual source identity attachment. Not a changed field, global cone, complete upstream norm gate, loop/repair/N or recursion.',
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Current original analytic Ra<R<=110 relaxed input and source attachment PASS',flush=True)
    return result


if __name__=='__main__':run()
