"""Check source-supported density formulas and original local integration.

Declared provider queries are not whole-chart or full-route admission.
The independent moderate fixtures exercise the unchanged density formulas.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_cutoff_density_oracle as current
import lei_ren_part1_paper_compliant_current_native_Rc_factored_source_oracle_check as inherited
import lei_ren_part1_paper_compliant_current_native_phase_first_jets_check as reference

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
cover=current.cover;density=current.density;packets=current.packets
require=reference.require;contains=reference.contains


def independent_density_checks(oracle):
    c=oracle.ctx;p=mp.mp.clone();p.dps=110;N=2048;count=0
    ledger=reference.qchecks.new_ledger();bases=tuple(c.mpf(0) for _ in range(5))
    def scalar(x):
        return current.cutoff.prior.ScaledEnclosure(current.cutoff.prior.FormalScale(bases),c.mpf(str(x)),ledger)
    for A0 in ('-.1','0','.1'):
        E,E_Z,V,V_Z=[scalar(x) for x in ('1.3','.3','.7','-.2')]
        primitives={key:scalar(x) for key,x in dict(A=A0,A_Z='.07',B_over_Pstar='.17',B_Z_over_Pstar='-.09').items()}
        old=density.density_Z_kernels(E,E_Z,V,V_Z,primitives,N)
        new=oracle.owner.kernel(E,E_Z,V,V_Z,primitives,N)
        e,ez,v,vz=[p.mpf(x) for x in ('1.3','.3','.7','-.2')]
        a=p.mpf(A0)/N;az=p.mpf('.07')/N;b=p.mpf('.17')/N;bz=p.mpf('-.09')/N
        de=e*p.expm1(a);dez=ez*p.expm1(a)+e*p.exp(a)*az
        expected=dict(kernels=dict(m=b,h=de,k=v*de+e*b+de*b,
            e=2*v*b+b*b-e*de-de*de/2,p=e*de+de*de/2),
            Z_derivatives=dict(m=bz,h=dez,k=vz*de+v*dez+ez*b+e*bz+dez*b+de*bz,
                e=2*vz*b+2*v*bz+2*b*bz-ez*de-e*dez-de*dez,
                p=ez*de+e*dez+de*dez))
        for group,rows in expected.items():
            for key,value in rows.items():
                target=c.mpf(p.nstr(value,115))
                for result in (old,new):
                    require(contains(cover.phase.bounded_value(result[group][key])+c.mpf(('-1e-95','1e-95')),target),
                        'Independent original density reference outside '+group+' '+key+' A='+A0)
                count+=1
        for name,value in (('original_E',E),('original_E_Z',E_Z),('original_V',V),('original_V_Z',V_Z)):
            require(new['velocities'][name] is value,'Supported factor adapter changed original E/V rows')
    return dict(passed=True,independent_signed_density_C0_Z_comparisons=count,
        signed_A_cases=['negative','zero','positive'],reference_precision_dps=110,
        original_and_supported_kernels_each_contain_independent_scalar_reference=True,
        synthetic_fixtures_not_original_field_values=True)


@current.native.inlet.source_precision
def run(role_owner,report=None):
    began=time.monotonic();built=role_owner.build();oracle=current.NativeCutoffFactoredOracle(role_owner,built)
    report=json.loads((HERE/current.NAME).read_bytes()) if report is None else report
    require(report[current.GATE] and report['source_family']==oracle.source_family,'Same original all17 producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed original source prerequisite: '+name)
    records=report['actual_original_declared_provider_queries']
    require(set(records)==set(current.original.POINTS) and len(records)==17,'All17 declared original providers required')
    require(report['enclosed_original_declared_queries']==17 and report['unresolved_original_declared_queries']==0,
        'All17 declared original provider queries must resolve')
    require(report['actual_original_role_dispatches']==160 and report['candidate_N']==2048 and report['full_Z_interval']==[-1,1],
        'Original declared N/Z/dispatch convention differs')
    for chart,row in records.items():
        require(row['status']=='enclosed','Original provider unresolved: '+chart)
        source=row['actual_original_spatial_source'];qsource=source['original_q_slow_jet_source']
        require(source['status']=='enclosed' and source['original_supported_density_AST_replacements']==[1,1],
            'Original supported density body adapter differs')
        require(source['nonlinear_density_precedes_cutoff_and_signed_u_branch_unions'] and
            source['overlapping_cutoff_or_signed_u_branches_hulled_not_added'],'Branch sum or premature nonlinear hull')
        require(source['original_A_over_N_formal_factor_and_A_Z_preserved'],'Original A/A_Z replaced by support bound')
        require(qsource['aggregate_source_roots_are_unconditioned_packet_metadata'] and
            qsource['nonlinear_phase_and_density_must_use_conditional_branch_query_objects'],'Aggregate source metadata not explicit')
        require(all(branch['spatial_query_status']=='enclosed' for branch in source['original_cutoff_branch_density_queries']),
            'A required cutoff branch was dropped')
        require(all(branch['original_source_query_AST_replacements']==[1,1]
            for branch in source['original_cutoff_branch_density_queries']),'Original conditional query adapter differs')
        require('separate_original_P0' in row and 'separate_original_P0_Z' in row,'Separate pressure datum lost')
    require(not report['all17_whole_chart_or_continuous_route_ranges_admitted'] and not any(report.get(k) for k in packets.OPEN),
        'Interior queries cannot close global gates')
    c=oracle.ctx;box=c.mpf((cover.ep(c.mpf('.12'))[0],cover.ep(c.mpf('.15'))[1]))
    broad=oracle.density_frame(chart='O2_slope',Z=(-1,1),coordinate=box,N=2048)
    require(broad.values is not None,'Original broad O2 source regressed')
    rows={(row.get('chart'),row.get('function_role')):row for row in built['graph'].nodes if row['operation']=='original_function_graph'}
    dispatches=0
    for key in density.RATES:
        for order in ('C0','Z'):
            got=oracle.dispatch_function_range(rows['O2_slope','density_'+key+'_'+order],broad)
            require(got is broad.values[order][key] and got.ctx is c,'Original range/context dispatch changed')
            dispatches+=1
    integral=oracle.owner.contribution(Z=(-1,1),left='.12',right='.15',N=2048)
    require(integral['record']['phase_and_overlapping_u_unions_before_single_integral_mass'],
        'Overlapping cutoff/phase cover must precede one positive integral mass')
    require(any(not v.zero for v in integral['contributions'].values()),'Original nonzero local source collapsed to zero')
    require(all(v.ctx is c for group in ('contributions','Z_derivatives') for v in integral[group].values()),
        'Original integral context changed')
    collar=oracle.density_frame(chart='bridge_first',Z=(-1,1),coordinate={'selected_sc_multiple':'1/2'},N=2048)
    require(collar.values is not None,'Original generic collar range remains unresolved')
    require('separate_original_P0' in collar.record and 'separate_original_P0_Z' in collar.record,'Collar pressure datum lost')
    require(collar.record['derived_original_global_phase']['phase_independent_of_Z'],'Original collar phase gained Z dependence')
    amplitude=oracle.amplitude_frame(Z=(-1,1),N=2048)
    for role in ('Rc_E_C0','Rc_E_Z'):
        got=oracle.dispatch_function_range(rows['O3_power',role],amplitude)
        require(got.ctx is c,'Separate Rc amplitude dispatch context changed');dispatches+=1
    independent=independent_density_checks(oracle)
    result=dict(all_passed=True,source_family=oracle.source_family,**{current.GATE:True},
        declared_original_provider_queries=17,enclosed_declared_original_provider_queries=17,
        unresolved_declared_original_provider_queries=0,declared_original_role_dispatches=160,
        independent_original_density_formula_checks=independent,
        inherited_role_namespace_root_checks=inherited.role_checks(role_owner,built,oracle),
        inherited_fail_closed_dispatch_checks=inherited.negative_checks(oracle,built,broad),
        additional_original_broad_cell_and_amplitude_dispatches=dispatches,
        actual_original_broad_O2_source=broad.record,actual_original_broad_O2_local_C0_Z_integral=integral['record'],
        actual_local_signed_integral_rows=10,original_generic_initial_collar_frame=collar.record,
        generic_collar_range_not_declared_exact_flat_or_zero=True,
        original_A_A_Z_formal_factors_and_nonzero_E_V_cross_terms_preserved=True,
        all17_whole_chart_or_continuous_route_ranges_admitted=False,full_factored_function_graph_evaluator_installed=False,
        actual_original_full_route_numerical_integrals_evaluated=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**oracle.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='All17 declared interior original provider ranges, independent density formulas, original broad O2 local C0/Z integration, generic collar range and separate Rc amplitude. Whole-route histories/targets, exact repair controls, global N, terminal closure and recursion remain open.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original cutoff density focused check PASS:17 providers,30 independent density comparisons',flush=True)
    return result
