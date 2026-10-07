"""Independent original paired Poisson derivatives and live source route."""
import ast
import inspect
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_paired_C1_transport as current
import lei_ren_part1_paper_compliant_current_native_periodic_C1_support_transport_check as accepted_checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
reference=accepted_checks.reference;require=reference.require


@current.native.inlet.source_precision
def independent_paired_parameter_checks(c):
    scales=reference.original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1e10',p2_abs_max='1e8',dps=120)
    p=scales.ctx;seed=reference.original.GenericShearLoop(scales,a='.8',b='.2',p1='1e9',p2=0,Utheta='1.3')
    def original(u,psi):
        loop=reference.original.GenericShearLoop(scales,a='.8',b='.2',p1='1e9',p2=u*scales.d_star/seed.q,Utheta='1.3')
        W1,W2=loop._w_integrals(psi)
        return W1/loop.h,W2/(loop.h*loop.h),loop
    count=0;cases=0;tol=p.mpf('1e-85')
    for utext in ('-1000','-1','-.2','0','.2','1','1000'):
        u=p.mpf(utext)
        for fraction in ('.001','.137','.49','.91','1'):
            psi=2*p.pi*p.mpf(fraction);P,H,loop=original(u,psi)
            Pu=p.diff(lambda v:original(v,psi)[0],u);Hu=p.diff(lambda v:original(v,psi)[1],u)
            r=loop.r;s=loop.one_minus_r2;D=loop._denominator(psi)
            if fraction=='1':E=2*p.pi
            else:
                minus=loop.one_minus_abs_r if r>=0 else 2-loop.one_minus_abs_r
                plus=2-loop.one_minus_abs_r if r>=0 else loop.one_minus_abs_r
                E=2*p.atan2(plus*p.sin(psi/2),minus*p.cos(psi/2))
            pairedP=P+u*Pu;pairedH=2*H+u*Hu
            expectedP=s**p.mpf('1.5')*p.sin(psi)/D
            expectedH=E+r*s*p.sin(psi)/D+s*s*p.sin(psi)*(p.cos(psi)-r)/(D*D)
            require(abs(pairedP-expectedP)<tol,'Independent original paired P identity failed')
            require(abs(pairedP)<=1+tol,'Independent original paired P bound failed')
            require(abs(pairedH-expectedH)<tol,'Independent original paired H identity failed')
            require(abs(pairedH)<=2*p.pi+3+tol,'Independent original paired H bound failed')
            require(abs(Pu)<=20*p.pi+tol,'Independent original P_u bound failed')
            require(abs(Hu)<=1000*p.pi+tol,'Independent original H_u bound failed')
            count+=6;cases+=1
    return dict(passed=True,independent_original_P_H_paired_identity_and_bound_comparisons=count,
        original_signed_small_zero_large_u_cases=cases,original_u_values=['-1000','-1','-.2','0','.2','1','1000'],
        original_GenericShearLoop_integrals_and_mpmath_diff_used=True,
        scalar_reference_precision_dps=120,comparison_allowance='1e-85',
        original_r_zero_uses_original_small_Fourier_limit=True,synthetic_fixtures_not_original_field_data=True)


@current.native.inlet.source_precision
def independent_paired_phase_checks(c):
    tree=ast.parse(inspect.getsource(accepted_checks.independent_periodic_checks))
    tree,counts=current.density.replace_expressions(tree,[
        ("current.periodic_C1_support(source['roots'],qrows,rows,c.ln(c.mpf('.7')),dstar)",
         "current.paired_C1_support(source['roots'],qrows,rows,c.ln(c.mpf('.7')),dstar)")])
    scope=dict(vars(accepted_checks));scope['current']=current
    exec(compile(tree,'<independent-original-scalar-reference-checks-of-paired-Z-support>','exec'),scope)
    got=scope['independent_periodic_checks'](c);got['paired_support_test_consumer_AST_replacements']=counts
    return got


@current.native.inlet.source_precision
def run(report,owner,live):
    began=time.monotonic()
    require(report[current.GATE] and report['source_family']==owner.family,'Same original paired Poisson producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed paired support prerequisite: '+name)
    checked=json.loads((HERE/current.accepted.RECEIPT).read_bytes())
    require(checked['all_passed'] and checked[current.accepted.GATE] and checked['source_family']==owner.family,
        'Accepted same-original whole-period support baseline required')
    require(report['enclosed_original_whole_cells']==24 and report['full24_original_C1_integral_range_transport_enclosed']
        and len(live['cells'])==24 and live['history'] is not None,'All24 original continuous first-Z transport cells required')
    integral_rows=0;decisions=0;used=0;retained=0;paired_traces=0
    for cell in live['cells']:
        for group in ('values','Z_derivatives'):
            for value in cell[group].values():
                require(value.ctx is owner.ctx and value.ledger is owner.coordinates.ledger,'Original integral arithmetic changed');integral_rows+=1
        if cell['frame'] is None:continue
        source=cell['frame'].record['actual_original_spatial_source']
        require(source['original_paired_Poisson_support_before_nonlinear_density'],'Original nonlinear ordering changed')
        for trace in source['original_whole_period_Z_support_ranges_before_nonlinear_density']:
            record=trace['original_whole_period_derivative_support'];theorem=record['original_paired_Poisson_parameter_theorem']
            require(theorem['passed'] and theorem['original_exact_paired_derivative_identities']==3 and
                record['original_whole_period_Z_support_AST_replacements']==[1,1,1] and
                record['original_linear_q_Z_and_nonzero_p2_Z_retained'],'Original paired derivative identity/source changed')
            require(trace['original_C0_A_B_objects_retained'] and trace['original_linear_q_Z_and_all_source_cross_rows_retained'],
                'Original C0 or source cross rows changed');paired_traces+=1
            for key in ('A_Z','B_Z_over_Pstar'):
                decisions+=1
                if trace['ordinary_Z_range_decisions'][key]['original_tighter_formal_range_retained']:retained+=1
                else:used+=1
    require(used>0 and integral_rows==240,'Actual source derivative supports and240 integral rows required')
    baseline=json.loads((HERE/current.accepted.NAME).read_bytes())
    for group in ('original_Rc_target_C0_ranges','original_Rc_target_Z_ranges'):
        for key,value in report[group].items():
            old=baseline[group][key]
            if old['exact_zero']:require(value['exact_zero'],'Original exact target zero changed')
            elif not value['exact_zero']:
                require(ep(value['log_absolute_upper'])[1]<=ep(packets.interval(owner.ctx,old['log_absolute_upper']))[1],
                    'Original target upper range worsened: '+group+' '+key)
    for row in report['original_serial_cells']:
        if row['chart']=='O3_power':
            require(row['incoming_C0']['p']==row['outgoing_C0']['p'] and row['incoming_Z']['p']==row['outgoing_Z']['p'],
                'Original quiet pressure memory changed')
    frontier=report['current_weighted_target_range_frontier'];oldfrontier=report['checked_unpaired_weighted_target_range_frontier']
    labels={cell['record']['label'] for cell in live['cells']}
    for data in (frontier,oldfrontier):
        require(data is not None and data['original_continuous_cells']==24 and len(data['weighted_original_cell_ranges'])==24,
            'Current and baseline actual weighted24-cell diagnostics required')
        require(data['applied_C0_and_Z_support_ranges_used'] and data['actual_density_primitive_ranges_inspected']==276,
            'Diagnostic must use actual applied density primitive supports')
        for key,dominant in data['dominant_normalized_target_contribution_cells'].items():
            require(dominant['label'] in labels,'Dominant contribution must refer to an actual original cell')
            maxima=[ep(cell['normalized_final_target_contribution_ranges'][key]['log_absolute_upper'])[1]
                for cell in data['weighted_original_cell_ranges'] if not cell['normalized_final_target_contribution_ranges'][key]['exact_zero']]
            require(ep(dominant['log_absolute_upper'])[1]==max(maxima),'Weighted dominant range rank incorrect')
    require(not any(report.get(key) for key in packets.OPEN) and not report['useful_repair_contraction_or_actual_controls_established'],
        'Paired parameter bounds cannot complete controls/global N/closure/recursion')
    independent=independent_paired_parameter_checks(owner.ctx);phase=independent_paired_phase_checks(owner.ctx)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=2048,
        original_continuous_route_cells=24,original_live_integral_C0_Z_rows=integral_rows,
        original_paired_first_jet_support_traces=paired_traces,original_branchwise_Z_range_decisions=decisions,
        tighter_original_paired_Z_support_ranges_used=used,original_tighter_Z_ranges_retained=retained,
        strict_target_absolute_upper_reductions=report['strict_target_absolute_upper_reductions'],
        independent_original_paired_parameter_checks=independent,independent_original_paired_phase_checks=phase,
        current_and_checked_baseline_weighted_source_range_frontiers_checked=True,
        original_linear_q_Z_and_C0_cross_terms_quiet_pressure_memory_retained=True,
        completed_actual_controls_global_N_closure_or_recursion=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Exact original paired Poisson bounds, independent scalar integral/first derivative references, actual24-cell transport and applied-support weighted target frontier. Genuine linear q_Z, actual controls/global N/closure/recursion remain open.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original paired Poisson focused check PASS;',used,'tighter derivative ranges;',report['strict_target_absolute_upper_reductions'],'strict target reductions',flush=True)
    return result
