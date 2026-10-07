"""Independent all-N density references and source-backed Rc target checks."""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets as current
import lei_ren_part1_paper_compliant_current_native_Rc_C1_histories_check as previous_checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;packets=current.packets;ep=current.ep
require=previous_checks.require;same=previous_checks.same;overlaps=previous_checks.overlaps


def independent_density_references(c):
    p=mp.mp.clone();p.dps=160;count=0;regions={};allow=c.mpf(('-1e-120','1e-120'))
    for label,A0,AZ,B0,BZ in (('signed_primitive_crossing','0','2','.5','-.1'),
        ('near_original_universal_cap','158','.5','-2','.25'),('tiny_nonzero_primitive','1e-40','2e-41','-5e-41','1e-41')):
        bases=tuple(c.mpf(0) for _ in range(5));ledger={}
        scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
        Z=c.mpf((-1,1));E=scalar(c.mpf('1.3')+c.mpf('.1')*Z);V=scalar(c.mpf('-.4')+c.mpf('.2')*Z)
        primitives=dict(A=scalar(c.mpf(A0)+c.mpf(AZ)*Z),A_Z=scalar(AZ),
            B_over_Pstar=scalar(c.mpf(B0)+c.mpf(BZ)*Z),B_Z_over_Pstar=scalar(BZ))
        got=current.uniform_density_orders(E,scalar('.1'),V,scalar('.2'),primitives)
        local=0
        for N in (160,1024,10**50):
            value=current.evaluate(got['values'],N);jet=current.evaluate(got['Z_derivatives'],N)
            def ref(z):
                e=p.mpf('1.3')+p.mpf('.1')*z;v=p.mpf('-.4')+p.mpf('.2')*z
                a=p.mpf(A0)+p.mpf(AZ)*z;b=p.mpf(B0)+p.mpf(BZ)*z
                changed_e=e*p.exp(a/N);changed_v=v+b/N
                # Independently subtract the original normalized five
                # history densities, rather than replaying the lift code.
                original=dict(m=v,h=e,k=e*v,e=v*v-e*e/2,p=e*e/2)
                changed=dict(m=changed_v,h=changed_e,k=changed_e*changed_v,
                    e=changed_v*changed_v-changed_e*changed_e/2,p=changed_e*changed_e/2)
                return {key:changed[key]-original[key] for key in current.RATES}
            for ztext in ('-.75','0','.75'):
                z=p.mpf(ztext);reference=ref(z)
                for key in current.RATES:
                    derivative=p.diff(lambda zz:ref(zz)[key],z)
                    require(previous_checks.previous_checks.contains(value[key].finite_interval()+allow,c.mpf(p.nstr(reference[key],165))),
                        'Independent all-N original density value not enclosed: '+label+' '+key)
                    require(previous_checks.previous_checks.contains(jet[key].finite_interval()+allow,c.mpf(p.nstr(derivative,165))),
                        'Independent all-N original density Z not enclosed: '+label+' '+key)
                    count+=2;local+=2
            if label=='tiny_nonzero_primitive':
                require(not value['h'].zero and not jet['h'].zero,'All-N lift lost tiny nonzero original expm1')
        regions[label]=dict(independent_original_density_C0_Z_comparisons=local,
            integer_N_cases=[160,1024,'10^50'],references_do_not_define_native_field_values=True)
    return dict(passed=True,independent_original_signed_density_C0_Z_comparisons=count,regions=regions,
        tiny_expm1_order_retained=True,signed_primitive_and_near_cap_cases_checked=True)


def record_value(row,coordinates):
    c=coordinates.ctx;s=row['formal_positive_scale'];powers=tuple(s['source_exponents'])+(s['radius_power'],)
    return current.prior.ScaledEnclosure(current.prior.FormalScale(coordinates.bases,powers,
        packets.interval(c,s['additional_log_interval'])),packets.interval(c,row['coefficient_interval']),coordinates.ledger)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2 and saved['all_integer_N_lower']==160,
        'Native source-backed all-N>=160 C0/Z target family required')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed original parameter target source: '+name)
    require(not saved['actual_control_functions_or_terminal_closure_installed']
        and not saved['one_global_finite_N_or_control_field_admitted'] and not any(saved.get(k) for k in packets.OPEN),
        'All-N source orders and repair conditions cannot install controls/global admission')
    old_manifest=json.loads((HERE/current.preceding.NAME).read_bytes())
    receipt=json.loads((HERE/current.preceding.RECEIPT).read_bytes())
    require(receipt['all_passed'] and receipt[current.preceding.GATE],'Original no-gap Rc source receipt required')
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeRcParameterTargets(current.preceding.NativeRcC1Histories(
            current.preceding.preceding.NativeO2C1Histories(current.preceding.preceding.make_middle_owner(bridge))))
        independent=independent_density_references(owner.ctx)
        require(saved['exact_density_order_theorem']==current.exact_density_order_theorem(),'Exact all-N signed density theorem changed')
        rows=0;inherited_rows=0;target_count=0;quiet_count=0;pressure_count=0;comparisons=0;regions={};bound_comparison={}
        expected=[dict(label=label,chart=chart,left=left,right=right) for label,chart,left,right in current.ROUTE]
        for name,record in saved['actual_original_Rc_parameter_target_records'].items():
            got=owner.route(packets.interval(owner.ctx,record['Z_box']))
            require(packets.encode(got['record'])==record,'Native all-N Rc target route changed: '+name)
            require(record['actual_original_route']==expected and record['original_true_cell_count']==24
                and record['complete_original_chart_count']==17 and len(record['original_seams'])==16,
                'All24 original cells/17charts/16seams required')
            incoming={key:{power:owner.coordinates.scalar(0) for power in current.ORDERS} for key in current.RATES}
            incomingZ={key:dict(row) for key,row in incoming.items()}
            for cell in got['cells']:
                label=cell['record']['label'];chart=cell['record']['chart'];geometry=cell['geometry']
                require(geometry['record']['width_and_endpoints_independent_of_Z'] and ep(geometry['width'].coefficient)[0]>0,
                    'Original positive fixed-Z log width required: '+label)
                if label=='initial_flat_collar':
                    require(cell['record']['original_signed_source_C1']['original_cutoff_and_primitive_C0_Z_exact_zero'],
                        'Original whole left collar exact-zero support required')
                else:
                    proof=cell['record']['original_signed_source_C1']['original_new_periodic_C1_cover']
                    require(proof['original_periodic_parameter_C1_theorem']==current.serial.periodic_parameter_theorem(),
                        'New parameter-uniform bounds must apply to every original active cell')
                    lift=cell['density']['record']
                    require(lift['uniform_N_lower']==160 and lift['coefficients_remain_N_dependent_functions']
                        and lift['original_universal_primitive_cap159_intersected_as_function_theorem'],
                        'All-N coefficient covers cannot become N-independent functions or clipped fields')
                    argument=packets.interval(owner.ctx,lift['primitive_over_N_union']);universal=owner.ctx.mpf((-159,159))/160
                    require(ep(argument)[0]>=ep(universal)[0] and ep(argument)[1]<=ep(universal)[1],
                        'Original universal source exp/g bound must be retained')
                for key in current.RATES:
                    decay=cell['factors'][key]['decay']
                    for power in current.ORDERS:
                        same(cell['incoming'][key][power],incoming[key][power],'Actual all-N C0 memory lost: '+label)
                        same(cell['incoming_Z'][key][power],incomingZ[key][power],'Actual all-N Z memory lost: '+label)
                        same(cell['cumulative'][key][power],decay*incoming[key][power]+cell['values'][key][power],
                            'Original all-N signed C0 transport mismatch: '+label)
                        same(cell['cumulative_Z'][key][power],decay*incomingZ[key][power]+cell['Z_derivatives'][key][power],
                            'Original all-N signed Z transport mismatch: '+label)
                        if chart=='O3_power':
                            require(cell['values'][key][power].zero and cell['Z_derivatives'][key][power].zero,
                                'Whole original quiet power must keep every N-order increment exactly zero');quiet_count+=2
                        for v in (cell['values'][key][power],cell['Z_derivatives'][key][power],cell['cumulative'][key][power],cell['cumulative_Z'][key][power]):
                            require(v.scale.bases is owner.coordinates.bases and v.ledger is owner.coordinates.ledger,
                                'All-N coefficients must keep original common bases/ledger')
                        rows+=4;inherited_rows+=2
                    if chart=='O3_power' and key=='p':
                        for power in current.ORDERS:
                            same(cell['cumulative'][key][power],incoming[key][power],'Quiet all-N pressure C0 memory lost')
                            same(cell['cumulative_Z'][key][power],incomingZ[key][power],'Quiet all-N pressure Z memory lost');pressure_count+=2
                incoming,incomingZ=cell['cumulative'],cell['cumulative_Z']
            amplitude=record['original_Rc_amplitude'];provenance=amplitude['same_original_source_provenance']
            require(provenance['chart']=='O3_power' and amplitude['checked_native_and_exact_repair_original_mu_equal'],
                'Actual Rc amplitude and exact repair must share original mu/source')
            require(packets.interval(owner.ctx,provenance['coordinate_box'])._mpi_==
                (2/owner.transfer.geometry.binder.fixed['Tw'])._mpi_,'Same original Rc E/E_Z source atphase2/Tw required')
            target=current.target_rows(got['values'],got['Z_derivatives'],got['amplitude'],got['amplitude_Z'],
                got['logA'],got['mu'],got['logmu'])
            for key in current.repair.ROWS:
                for power in current.ORDERS:
                    same(target['values'][key][power],got['targets']['values'][key][power],'Original signed target coefficient mismatch')
                    same(target['Z_derivatives'][key][power],got['targets']['Z_derivatives'][key][power],'Original signed target Z mismatch');target_count+=2
            contract=record['original_joint_target_function_contract']
            require(contract['joint_divided_row_quantitative_cancellation_not_established']
                and contract['common_basis_does_not_restore_lost_interval_correlations']
                and contract['exact_original_implicit_equation']=='B_exact(mu)*h+N*r(N,Z)+Q_exact(mu,h)/N=0'
                and contract['P0_and_P0_Z_unchanged'] and not contract['actual_control_functions_installed'],
                'Original divided row/sign/P0 and uninstalled controls must remain explicit')
            require(packets.encode(current.target_caps(target))==record['actual_uniform_N_scaled_repair_C1_caps'],
                'Actual native all-N target caps changed')
            conditions=got['conditions']
            require(conditions['contraction_at_most']=='1/2' and conditions['image_radius_at_most']=='3*rho/4',
                'Checked exact generic inverse condition recipe required')
            require(ep(conditions['repair_sufficient_common_log_N_lower'])[0]>=ep(owner.ctx.ln(160))[0]
                and record['no_log_N_exponential_or_huge_integer_materialized'],
                'Native repair log condition must include original all-N floor without materializing N')
            require(len(record['actual_weighted_cell_target_diagnostics'])==24
                and set(record['dominant_target_bound_cell_labels'])==set(current.repair.ROWS),
                'Native source target budgets must identify their dominating cells')
            new=current.evaluate(got['values'],1024);newZ=current.evaluate(got['Z_derivatives'],1024)
            old=old_manifest['actual_original_inlet_to_Rc_C1_records'][name];diagnostic={}
            for key in current.RATES:
                oldv=record_value(old['actual_Rc_correction_C0'][key],owner.coordinates)
                oldz=record_value(old['actual_Rc_correction_Z'][key],owner.coordinates)
                require(overlaps(new[key],oldv) and overlaps(newZ[key],oldz),
                    'Same original N1024 source covers disagree after uniform lift: '+key);comparisons+=2
                diagnostic[key]=dict(C0_new_upper_log=current.magnitude(new[key]).record(),
                    C0_prior_upper_log=current.magnitude(oldv).record(),Z_new_upper_log=current.magnitude(newZ[key]).record(),
                    Z_prior_upper_log=current.magnitude(oldz).record(),bound_comparison_not_true_error_measurement=True)
            bound_comparison[name]=diagnostic
            regions[name]=dict(original_cells=24,original_charts=17,new_C0_Z_order_and_transport_rows=960,
                inherited_order_rows=480,normalized_repair_target_C0_Z_order_rows=20,
                whole_N_scaled_target_cap=record['actual_uniform_N_scaled_repair_C1_caps']['whole_target_C1_cap'],
                dominant_target_bound_cells=record['dominant_target_bound_cell_labels'])
            print('Actual all-N Rc targets checked:',name,flush=True)
        try:current.evaluate(got['values'],159)
        except ValueError:pass
        else:raise ArithmeticError('Unsafe all-N source frequency accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},native_Z_queries_checked=2,
        all_integer_N_lower=160,original_true_cells_checked=24,original_charts_checked=17,
        actual_C0_Z_order_and_transport_rows_checked=rows,actual_inherited_order_rows_checked=inherited_rows,
        actual_normalized_target_C0_Z_order_rows_checked=target_count,exact_quiet_power_N_order_rows_checked=quiet_count,
        preserved_quiet_pressure_N_order_rows_checked=pressure_count,
        same_original_N1024_Rc_cover_consistency_rows_checked=comparisons,
        independent_density_reference_checks=independent,exact_density_order_theorem=current.exact_density_order_theorem(),
        native_uniform_target_caps_connected_to_checked_exact_repair_inverse=True,
        joint_divided_row_and_same_original_mu_P0_provenance_checked=True,
        uniform_new_periodic_bounds_on_entire_actual_original_route_checked=True,
        same_N_prior_and_new_Rc_bound_diagnostics=bound_comparison,regions=regions,
        actual_control_functions_or_terminal_closure_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('All-N Rc target PASS:',rows,'order/transport rows;',target_count,'target rows',flush=True)
    return result


if __name__=='__main__':run()
