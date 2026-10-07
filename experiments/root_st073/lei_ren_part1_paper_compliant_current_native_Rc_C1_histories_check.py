"""Actual original Rc endpoint, quiet power and no-gap C0/Z memory checks."""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_C1_histories as current
import lei_ren_part1_paper_compliant_current_native_O2_C1_histories_check as previous_checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;packets=current.packets;ep=current.ep
require=previous_checks.require;same=previous_checks.same;overlaps=previous_checks.overlaps


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2 and saved['actual_added_O3_Rc_cell_count']==3,
        'Actual original three-cell O3/Rc history route required')
    require(saved['original_entire_inlet_to_Rc_radial_route_covered']
        and not saved['global_inlet_to_Rc_histories_admitted'] and not any(saved.get(k) for k in packets.OPEN),
        'No-gap conservative Rc covers cannot imply quantitative terminal/global admission')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed original Rc source: '+name)
    previous=json.loads((HERE/current.preceding.NAME).read_bytes())
    inherited=json.loads((HERE/current.preceding.RECEIPT).read_bytes())
    require(inherited['all_passed'] and inherited[current.preceding.GATE],'Accepted actual O3-inlet history source required')
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeRcC1Histories(current.preceding.NativeO2C1Histories(current.preceding.make_middle_owner(bridge)))
        require(saved['original_O3_transition_power_source_join']==owner.join,'Original transition/power source join changed')
        require(saved['original_Rc_reservation']==owner.reservation,'Original Rc reservation changed')
        Tw=owner.transfer.geometry.binder.fixed['Tw'];require(ep(Tw-2)[0]>0,'Actual Rc must lie inside original power chart')
        root_owner=owner.q_owner.owner.owner
        eta=packets.interval(owner.ctx,root_owner.scales['selected_positive_eta_log'])
        domain_mu=packets.interval(owner.ctx,owner.domain['right_collar_mu'])
        require(ep(eta)[1]<=ep(owner.ctx.ln(domain_mu))[0],
            'Original selected eta<=mu ensures entire canonical power is flat')
        rows=0;incoming_rows=0;quiet_rows=0;pressure_rows=0;composite_rows=0;terminal_rows=0;regions={}
        expected=[dict(label=label,chart=chart,left=left,right=right,
            coordinate_kind='original_power_offset' if chart=='O3_power' else 'original_transition_offset') for label,chart,left,right in current.ROUTE]
        for name,old in saved['actual_original_inlet_to_Rc_C1_records'].items():
            got=owner.route(packets.interval(owner.ctx,old['Z_box']),saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual original Rc C1 route changed: '+name)
            require(old['known_actual_original_O3_inlet_C1_history']==previous['actual_original_inlet_to_O2_exit_C1_records'][name],
                'Actual O3 coordinate0 incoming differs from accepted O2 exit')
            require(old['ordered_original_O3_Rc_cells']==expected and old['no_original_interval_skipped_from_true_inlet_to_actual_Rc'],
                'Original transition/power route must have no gap to actual Rc')
            incoming=got['initial']['correction']
            require(any(not v.zero for v in incoming['values'].values()) and any(not v.zero for v in incoming['Z_derivatives'].values()),
                'Actual inherited O3-inlet boxes cannot be zeroed')
            for index,(label,chart,left,right) in enumerate(current.ROUTE,1):
                cell=got['cells'][label];row=old['actual_O3_Rc_serial_C1_records'][label]
                geometry=cell['geometry']['record']
                require(geometry['width_and_endpoints_independent_of_Z'] and ep(cell['geometry']['width'].coefficient)==(1,1)
                    and ep(cell['geometry']['width'].scale.evaluate())==(0,0),
                    'Every added cell has exact original unit log-radius length')
                require(row['original_chart_uniform_a_positive_certificate']['whole_actual_source_positive_not_inferred_from_saved_denominator_box'],
                    'Original O3 source-function positivity theorem required')
                require(row['original_true_width_and_native_ordinary_y_conversion_applied_once']
                    and row['actual_incoming_correction_not_reset_or_background_double_added'],
                    'Original fixed-Z density/width/recovery path must be retained')
                require(row['O3_inlet_to_current_endpoint_C1_operator']['steps']==index,'Cumulative O3 operator lost a cell')
                if chart=='O3_slope_mu':
                    require(cell['geometry']['coordinate']._mpi_==owner.ctx.mpf((0,1))._mpi_,
                        'Whole original transition domain required')
                    require(row['original_chart_uniform_a_positive_certificate']['exact_a_minus2_source']=='2*mu*sigma(t)',
                        'Variable original transition shear cannot be fixed to power slope')
                else:
                    a=owner.ctx.mpf(left)/Tw;b=owner.ctx.mpf(right)/Tw
                    require(cell['geometry']['coordinate']._mpi_==owner.ctx.mpf((ep(a)[0],ep(b)[1]))._mpi_,
                        'Original full power-offset cell must use phase=offset/Tw')
                    require(row['original_right_endpoint']['original_power_offset']==right,
                        'Original power endpoint offset changed')
                    proof=row['original_periodic_parameter_C1_cover']
                    require(proof['original_full_box_q_and_q_Z_exact_zero_implies_primitive_C0_Z_exact_zero'],
                        'Entire original power source must have checked flat support')
                    roots=cell['source']['source']['roots'];positive=row['original_chart_uniform_a_positive_certificate']
                    mu=packets.interval(owner.ctx,positive['actual_positive_mu'])
                    require(mu._mpi_==domain_mu._mpi_ and positive['exact_a_minus2_source']=='2*mu',
                        'Same original positive canonical mu source required')
                    require(roots['b'][current.ZERO].zero and roots['b'][current.DZ].zero
                        and roots['a'][current.DZ].zero and roots['kappa_minus2'][current.DZ].zero,
                        'Original b=0 and Z-independent canonical power shear must be retained')
                    same(roots['kappa_minus2'][current.ZERO],roots['kappa_minus2'][current.ZERO].scalar(2*mu),
                        'Canonical tiny 2mu must not be lost through rounded a-2 subtraction')
                    require(all(v.zero for v in cell['primitives']['values'].values()),'Quiet power primitive first-Z must be exact zero')
                    for key in current.RATES:
                        require(cell['kernels']['kernels'][key].zero and cell['kernels']['Z_derivatives'][key].zero
                            and cell['contributions'][key].zero and cell['Z_derivatives'][key].zero,
                            'Quiet power density and true-width increments must be zero: '+key);quiet_rows+=2
                    endpoint=row['original_right_background_and_separate_P0_Z']['source_provenance']
                    require(endpoint['chart']=='O3_power' and packets.interval(owner.ctx,endpoint['coordinate_box'])._mpi_==b._mpi_,
                        'Original endpoint background must be directly queried at its offset/Tw')
                if label=='power_to_r_plus':
                    require(row['original_background_history_P0_join_receipt']==owner.join
                        and row['original_same_radius_phase_or_same_chart_seam']=='O3_slope_mu -> O3_power',
                        'Original Rw five-history/P0 source join must be bound')
                    require(row['original_same_radius_phase_or_same_chart_seam'] in
                        geometry['exact_original_radius_Jacobian_identities']['original_same_radius_periodic_phase_seam_identities'],
                        'Exact original transition/power radius/phase seam required')
                if label=='power_to_Rc':
                    require(row['original_background_history_P0_join_receipt']['same_original_chart_and_source_function'],
                        'Same original power source required across r_plus')
                for key in current.RATES:
                    same(cell['incoming']['values'][key],incoming['values'][key],'Actual inherited C0 lost: '+label+' '+key)
                    same(cell['incoming']['Z_derivatives'][key],incoming['Z_derivatives'][key],'Actual inherited Z lost: '+label+' '+key)
                    decay=cell['factors'][key]['decay']
                    same(cell['correction']['values'][key],decay*incoming['values'][key]+cell['contributions'][key],
                        'Original signed C0 transfer mismatch: '+label+' '+key)
                    same(cell['correction']['Z_derivatives'][key],decay*incoming['Z_derivatives'][key]+cell['Z_derivatives'][key],
                        'Original signed Z transfer mismatch: '+label+' '+key)
                    same(cell['own'][key],cell['background']['originals'][key]+cell['correction']['values'][key],
                        'Original endpoint own C0 mismatch: '+label+' '+key)
                    same(cell['own_Z'][key],cell['background']['Z_derivatives'][key]+cell['correction']['Z_derivatives'][key],
                        'Original endpoint own Z mismatch: '+label+' '+key)
                require(ep(cell['factors']['p']['decay'].coefficient)==(1,1) and ep(cell['factors']['p']['decay'].scale.evaluate())==(0,0),
                    'Rate-zero pressure must retain coefficient1')
                if chart=='O3_power':
                    same(cell['correction']['values']['p'],incoming['values']['p'],'Quiet power lost actual pressure memory')
                    same(cell['correction']['Z_derivatives']['p'],incoming['Z_derivatives']['p'],'Quiet power lost actual pressure Z memory');pressure_rows+=2
                require(cell['background']['record']['P0_not_merged_into_pressure_history'],'Original P0 and P0_Z must stay separate')
                for value in [*cell['own'].values(),*cell['own_Z'].values(),*cell['correction']['values'].values(),*cell['correction']['Z_derivatives'].values()]:
                    require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,
                        'Rc functions must share the accepted common arithmetic bases/ledger')
                require(not row['global_inlet_to_Rc_histories_admitted'] and not any(row.get(k) for k in packets.OPEN),
                    'Conservative Rc cells cannot admit quantitative global gates')
                rows+=20;incoming_rows+=10;incoming=cell['correction']
            composite=got['cumulative'].apply(got['initial']['correction']['values'],got['initial']['correction']['Z_derivatives'],owner.family)
            for key in current.RATES:
                for kind in ('values','Z_derivatives'):
                    require(overlaps(composite[kind][key],got['correction'][kind][key]),
                        'Whole O3 composite and serial covers disagree: '+key);composite_rows+=1
            endpoint=old['original_Rc_endpoint']
            require(endpoint['radius_expression']=='Rw*exp(2)' and endpoint['original_power_offset']==2
                and endpoint['original_power_phase_expression']=='2/Tw' and endpoint['original_Rc_endpoint_not_whole_power_phase1'],
                'Actual Rc must remain original power offset2, not phase1')
            require(packets.interval(owner.ctx,endpoint['actual_direct_endpoint_phase_cover'])._mpi_==(2/Tw)._mpi_,
                'Actual Rc phase query does not use the original source Tw')
            collar=old['original_right_collar_inside_power']
            require(packets.interval(owner.ctx,collar['original_power_offset_box'])._mpi_==owner.ctx.mpf(('.75','1.25'))._mpi_
                and collar['whole_original_power_q_flat_source_proof_covers_both_sides'],
                'Original two-sided right collar must be inside the power coverage')
            for key in current.RATES:
                same(got['own'][key],got['background']['originals'][key]+got['correction']['values'][key],
                    'Actual Rc own C0 mismatch: '+key)
                same(got['own_Z'][key],got['background']['Z_derivatives'][key]+got['correction']['Z_derivatives'][key],
                    'Actual Rc own Z mismatch: '+key);terminal_rows+=2
            require(old['quantitative_terminal_five_Z_identities_not_established'] and old['original_outer_2Rc_to_Rb_admission_still_missing'],
                'Terminal repair and outer admission remain unresolved')
            regions[name]=dict(actual_O3_Rc_cells=3,actual_correction_and_own_C0_Z_rows=60,
                actual_inherited_C0_Z_rows=30,exact_quiet_power_density_C0_Z_rows=20,
                preserved_pressure_C0_Z_rows=4,serial_composite_overlap_rows=10,actual_Rc_own_C0_Z_rows=10)
            print('Actual original Rc C1 covers checked:',name,flush=True)
        try:owner.route(N=159)
        except ValueError:pass
        else:raise ArithmeticError('Unsafe original whole-period frequency accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        native_Z_queries_checked=2,actual_added_O3_Rc_cells_checked=3,
        actual_added_cell_correction_and_own_C0_Z_rows_checked=rows,actual_inherited_C0_Z_rows_checked=incoming_rows,
        original_power_exact_zero_density_C0_Z_rows_checked=quiet_rows,
        original_quiet_power_preserved_pressure_C0_Z_rows_checked=pressure_rows,
        O3_serial_composite_overlap_rows_checked=composite_rows,actual_Rc_own_C0_Z_rows_checked=terminal_rows,
        actual_original_Rc_power_offset2_phase2_over_Tw_direct_endpoint_checked=True,
        original_two_sided_right_collar_inside_quiet_power_checked=True,
        actual_original_inlet_to_Rc_no_gap_function_domain_C0_Z_covers_checked=True,
        accepted_original_periodic_parameter_C1_and_scalar_reference_checks_inherited=True,
        original_transition_power_source_five_history_and_separate_P0_join_bound=True,
        conservative_covers_not_quantitative_terminal_closure=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),regions=regions,
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rc C1 PASS:',rows,'actual cell rows;',quiet_rows,'exact quiet density rows',flush=True)
    return result


if __name__=='__main__':run()
