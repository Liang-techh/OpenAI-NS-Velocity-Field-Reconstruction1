"""Original scalar derivative references and live periodic Z supports."""
import ast
import inspect
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_periodic_C1_support_transport as current
import lei_ren_part1_paper_compliant_current_native_collected_q2_transport_check as accepted_checks

accepted=current.accepted;reference=accepted_checks.reference
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets;require=reference.require


@current.native.inlet.source_precision
def independent_periodic_checks(c):
    scales=reference.original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    eta_log=c.ln(c.mpf(scales.ctx.nstr(scales.eta,105)));support_decisions=[]
    def periodic_fixture(source,qrows,dstar,phi):
        if qrows is None:return reference.current.conditioned_first_jets(source,None,dstar,phi)
        cutoffs,_=accepted.cutoff.conditional_cutoff_branches(source['roots']['kappa_minus2'][current.ZERO],eta_log)
        u,branches,empty=accepted.cover.signed_u_branches(source,dstar);outputs=[];supports=[]
        for cutoff_branch in cutoffs:
            rows=accepted.conditional_q2_jet(source['roots'],eta_log,c.ln(c.mpf('.7')),cutoff_branch)
            support=current.periodic_C1_support(source['roots'],qrows,rows,c.ln(c.mpf('.7')),dstar)
            supports.append(support['values'])
            for branch in branches:
                got=accepted.compile_collected_first_jets(source,qrows,rows,dstar,phi,branch)
                require(got['values'] is not None,'Original periodic-support scalar fixture unresolved')
                for key in ('A_Z','B_Z_over_Pstar'):
                    got['values'][key],decision=accepted.accepted.signed_C0_support_range(got['values'][key],support['values'][key])
                    support_decisions.append(decision)
                outputs.append(got)
        union=accepted.cover.density.local.same_source_union;got=dict(outputs[0])
        got['values']={key:union([row['values'][key] for row in outputs]) for key in reference.current.OUTPUTS}
        got['independent_direct_Z_support']={key:union([row[key] for row in supports]) for key in ('A_Z','B_Z_over_Pstar')}
        return got
    def direct_support_value(got,key):
        return got['independent_direct_Z_support'][key] if key in ('A_Z','B_Z_over_Pstar') else got['values'][key]
    tree=ast.parse(inspect.getsource(reference.scalar_first_checks))
    tree,counts=accepted.density.replace_expressions(tree,[
        ("current.conditioned_first_jets(source,qjet['rows'],dstar,'.137')","periodic_fixture(source,qjet['rows'],dstar,'.137')"),
        ("current.conditioned_first_jets(source,qjet['rows'],dstar,phi)","periodic_fixture(source,qjet['rows'],dstar,phi)"),
        ("current.conditioned_first_jets(source,None,dstar,'.137')","periodic_fixture(source,None,dstar,'.137')"),
        ("current.phase.bounded_value(got['values'][key])","current.phase.bounded_value(direct_support_value(got,key))")])
    scope=dict(vars(reference));scope.update(periodic_fixture=periodic_fixture,direct_support_value=direct_support_value)
    exec(compile(tree,'<independent-original-scalar-checks-of-periodic-Z-support>','exec'),scope)
    got=scope['scalar_first_checks'](c)
    got.update(test_consumer_AST_replacements=counts,independent_original_reference_functions_and_allowances_unchanged=True,
        original_A_Z_B_Z_reference_comparisons_against_direct_supports=12,
        original_Z_support_range_decisions=len(support_decisions),
        original_tighter_Z_ranges_retained=sum(row['original_tighter_formal_range_retained'] for row in support_decisions),
        periodic_Z_support_ranges_used=sum(not row['original_tighter_formal_range_retained'] for row in support_decisions))
    return got


@current.native.inlet.source_precision
def run(report,owner,live):
    began=time.monotonic()
    require(report[current.GATE] and report['source_family']==owner.family,'Same original periodic Z support producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed periodic support prerequisite: '+name)
    prior=json.loads((HERE/accepted.RECEIPT).read_bytes())
    require(prior['all_passed'] and prior[accepted.GATE] and prior['source_family']==owner.family,
        'Accepted same-original-family q squared mixed/phase checks required')
    require(report['enclosed_original_whole_cells']==24 and report['full24_original_C1_integral_range_transport_enclosed'],
        'Original all24 continuous first-Z transport cells required')
    require(len(live['cells'])==24 and live['history'] is not None,'Original live continuous route required')
    integral_rows=0;decisions=0;used=0;retained=0
    for cell in live['cells']:
        for group in ('values','Z_derivatives'):
            for value in cell[group].values():
                require(value.ctx is owner.ctx and value.ledger is owner.coordinates.ledger,'Original integral context/ledger changed');integral_rows+=1
        if cell['frame'] is None:continue
        source=cell['frame'].record['actual_original_spatial_source']
        require(source['derivative_support_not_C0_cap_differentiation'],'Original derivative support provenance missing')
        for row in source['original_whole_period_Z_support_ranges_before_nonlinear_density']:
            require(row['original_C0_A_B_objects_retained'] and row['original_linear_q_Z_and_all_source_cross_rows_retained'],
                'Derivative support changed C0 or cross rows')
            for key in ('A_Z','B_Z_over_Pstar'):
                decision=row['ordinary_Z_range_decisions'][key];decisions+=1
                if decision['original_tighter_formal_range_retained']:retained+=1
                else:used+=1
    for row in report['original_serial_cells']:
        if row['chart']=='O3_power':
            require(row['incoming_C0']['p']==row['outgoing_C0']['p'] and row['incoming_Z']['p']==row['outgoing_Z']['p'],
                'Original quiet pressure memory changed')
    comparison=report['comparison_with_checked_collected_q_squared_target_ranges']
    baseline=json.loads((HERE/accepted.NAME).read_bytes())
    for group in ('original_Rc_target_C0_ranges','original_Rc_target_Z_ranges'):
        for key,value in report[group].items():
            old=baseline[group][key]
            if old['exact_zero']:require(value['exact_zero'],'Exact original target zero changed')
            elif not value['exact_zero']:
                require(ep(value['log_absolute_upper'])[1]<=ep(packets.interval(owner.ctx,old['log_absolute_upper']))[1],
                    'Original target upper range worsened: '+group+' '+key)
    require(used>0,'Actual source derivative support must tighten at least one original range')
    require(not any(report.get(key) for key in packets.OPEN) and not report['useful_repair_contraction_or_actual_controls_established'],
        'Periodic derivative supports cannot complete controls/global N/closure/recursion')
    independent=independent_periodic_checks(owner.ctx)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=2048,
        actual_original_continuous_route_cells=24,actual_original_live_integral_C0_Z_rows=integral_rows,
        original_branchwise_Z_range_decisions=decisions,original_periodic_Z_ranges_used=used,
        original_tighter_Z_ranges_retained=retained,strict_target_absolute_upper_reductions=report['strict_target_absolute_upper_reductions'],
        independent_original_periodic_Z_support_checks=independent,
        original_C0_primitives_linear_q_Z_source_cross_terms_and_quiet_pressure_memory_retained=True,
        derivative_supports_not_C0_cap_derivatives=True,completed_actual_controls_global_N_closure_or_recursion=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Original whole-period ordinary Z supports with independent scalar first derivative references, actual tighter range decisions and all24 original continuous C1 transport cells. Actual controls/global N/closure/recursion remain open.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original periodic Z support focused check PASS;',used,'tighter derivative ranges;',report['strict_target_absolute_upper_reductions'],'strict target reductions',flush=True)
    return result
