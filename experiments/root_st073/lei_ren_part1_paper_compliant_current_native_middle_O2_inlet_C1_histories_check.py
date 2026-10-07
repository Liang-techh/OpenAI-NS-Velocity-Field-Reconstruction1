"""Focused no-gap middle/O2-inlet C1 route and inherited memory verification.

Accepted original loop, mass and affine theorems are inherited. New native
full-box source queries and actual route inputs/outputs are checked here.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_middle_O2_inlet_C1_histories as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;packets=current.packets
require=checks.require;ep=current.ep


def same(a,b,message):require(packets.encode(a.record())==packets.encode(b.record()),message)


def overlaps(a,b):
    """Compare both formal signed covers in one directed arithmetic coordinate."""
    c=a.ctx;a.scale.pair(b.scale)
    ref=c.mpf(max(ep(a.scale.evaluate())[1],ep(b.scale.evaluate())[1]))
    aa=a.coefficient*a.bounded_exp(a.scale.evaluate()-ref)
    bb=b.coefficient*b.bounded_exp(b.scale.evaluate()-ref)
    return max(ep(aa)[0],ep(bb)[0])<=min(ep(aa)[1],ep(bb)[1])


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2 and saved['complete_added_original_chart_count']==6,
        'Actual six-chart whole middle/O2-inlet C1 route required')
    require(not any(saved.get(k) for k in packets.OPEN),'Middle/O2 inlet covers cannot admit global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed middle/O2-inlet source: '+name)
    previous=json.loads((HERE/current.preceding.NAME).read_bytes())
    inherited=json.loads((HERE/current.preceding.RECEIPT).read_bytes())
    require(inherited['all_passed'] and inherited[current.preceding.GATE],'Accepted actual R110 history receipt required')
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        p=current.p
        c1=current.density.NativeDensityC1LocalIntegrals(p.first.NativePhaseFirstJets(p.slow.NativeQSlowJets(
            p.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))))
        second=current.preceding.preceding.NativeSecondBridgeC1Histories(p.NativeFirstBridgeC1Histories(
            current.transfer.NativeTrueChartC1Transfer(current.history.NativeC1HistoryTransfer(c1))))
        owner=current.NativeMiddleO2InletC1Histories(current.preceding.NativeBridgeSwitchC1Histories(second))
        rows=0;incoming_rows=0;overlap_rows=0;terminal_rows=0;regions={}
        expected=[dict(chart=chart,left=left,right=right) for chart,left,right in current.ROUTE]
        for name,old in saved['actual_inlet_to_O2_inlet_C1_history_records'].items():
            got=owner.route(packets.interval(owner.ctx,old['Z_box']),saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual middle/O2-inlet C1 route changed: '+name)
            require(old['known_actual_inlet_to_R110_C1_history']==previous['actual_inlet_to_R110_C1_history_records'][name],
                'Actual R110 incoming changed from accepted source')
            require(old['ordered_original_native_route']==expected and old['no_original_interval_skipped_between_true_inlet_and_O2_inlet'],
                'Original complete ordered six-chart route required')
            incoming=got['initial']['correction'];previous_chart,previous_endpoint='switch_power',1
            require(any(not v.zero for v in incoming['values'].values()),'Phase2 correction boxes cannot be replaced by zero')
            require(any(not v.zero for v in incoming['Z_derivatives'].values()),'Phase2 first-Z boxes cannot be replaced by zero')
            for index,(chart,left,right) in enumerate(current.ROUTE,1):
                cell=got['cells'][chart];row=old['actual_serial_chart_C1_history_records'][chart]
                require(row['original_left_endpoint']==dict(chart=previous_chart,coordinate=previous_endpoint),
                    'Wrong inherited original chart endpoint: '+chart)
                if chart=='actual_patch':
                    analytic=owner.ctx.exp(1)
                    require(ep(cell['geometry']['coordinate'])==(1,ep(analytic)[1]),'Whole analytic patch domain required')
                    require(ep(cell['geometry']['width'].coefficient)==(1,1) and ep(cell['geometry']['width'].scale.evaluate())==(0,0),
                        'Full patch log-radius width must be exactly1')
                    require(row['original_right_endpoint']['analytic_source_expression']=='exp(1)',
                        'Patch endpoint must retain the original analytic expression')
                    endpoint=packets.interval(owner.ctx,row['actual_right_original_background_and_separate_P0_Z']['source_provenance']['coordinate_box'])
                    require(endpoint._mpi_==analytic._mpi_,'Patch background must be queried at exp(1), not the full source box')
                else:
                    require(ep(cell['geometry']['coordinate'])==(left,right),'Entire native coordinate domain required: '+chart)
                join=row['original_background_history_P0_source_join_receipt']
                require(join==owner.joins[chart] and join['background_history_and_analytic_P0_identity_not_inferred_from_radius_only'],
                    'Accepted same-function background/history/P0 join required: '+chart)
                if chart=='actual_patch':require(row['original_actual_patch_direct_source_receipt']==owner.joins['Rh_reference'],
                    'Direct accepted original patch source admission required')
                geometry=cell['geometry']['record']
                require(geometry['width_and_endpoints_independent_of_Z'] and not cell['geometry']['width'].zero
                    and ep(cell['geometry']['width'].coefficient)[0]>0,'Original positive fixed-Z width required: '+chart)
                seam=previous_chart+' -> '+chart
                require(row['exact_original_neighbor_radius_and_periodic_phase_seam']==seam
                    and seam in geometry['exact_original_radius_Jacobian_identities']['original_same_radius_periodic_phase_seam_identities'],
                    'Exact original radius/phase seam absent: '+chart)
                require(row['original_chart_uniform_a_positive_certificate']['source_function_positivity_not_inferred_from_saved_box'],
                    'Point-source box cannot define positive theorem: '+chart)
                require(row['original_source_rows_are_ordinary_log_radius_and_Z'] and row['original_true_width_and_source_Jacobian_applied_once'],
                    'Original ordinary-y normalization must be retained: '+chart)
                require(row['R110_to_current_endpoint_C1_operator']['steps']==index,'Cumulative affine operator lost a cell')
                for key in current.RATES:
                    same(cell['incoming']['values'][key],incoming['values'][key],'Actual previous C0 memory not inherited: '+chart+' '+key)
                    same(cell['incoming']['Z_derivatives'][key],incoming['Z_derivatives'][key],'Actual previous Z memory not inherited: '+chart+' '+key)
                    decay=cell['factors'][key]['decay']
                    same(cell['correction']['values'][key],decay*incoming['values'][key]+cell['contributions'][key],
                        'Signed nonzero-incoming C0 transfer mismatch: '+chart+' '+key)
                    same(cell['correction']['Z_derivatives'][key],decay*incoming['Z_derivatives'][key]+cell['Z_derivatives'][key],
                        'Signed incoming Z transfer mismatch: '+chart+' '+key)
                    same(cell['own'][key],cell['background']['originals'][key]+cell['correction']['values'][key],
                        'Original endpoint background reset or double-added: '+chart+' '+key)
                    same(cell['own_Z'][key],cell['background']['Z_derivatives'][key]+cell['correction']['Z_derivatives'][key],
                        'Original endpoint background Z reset or double-added: '+chart+' '+key)
                require(ep(cell['factors']['p']['decay'].coefficient)==(1,1) and ep(cell['factors']['p']['decay'].scale.evaluate())==(0,0),
                    'Rate-zero pressure memory must survive every chart exactly')
                require(row['actual_right_original_background_and_separate_P0_Z']['P0_not_merged_into_pressure_history'],
                    'Analytic P0/P0_Z must stay separate: '+chart)
                for value in [*cell['own'].values(),*cell['own_Z'].values(),*cell['correction']['values'].values(),*cell['correction']['Z_derivatives'].values()]:
                    require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,
                        'Every actual C0/Z function must use the same common basis/ledger')
                require(not row['global_inlet_to_Rc_histories_admitted'] and not any(row.get(k) for k in packets.OPEN),
                    'Local route cannot complete quantitative global gates')
                rows+=20;incoming_rows+=10;incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
            composite=got['cumulative'].apply(got['initial']['correction']['values'],got['initial']['correction']['Z_derivatives'],owner.family)
            for key in current.RATES:
                require(overlaps(composite['values'][key],got['correction']['values'][key]),'Cumulative and serial C0 covers disagree: '+key)
                require(overlaps(composite['Z_derivatives'][key],got['correction']['Z_derivatives'][key]),'Cumulative and serial Z covers disagree: '+key)
                overlap_rows+=2
            require(got['cumulative'].steps==6 and old['original_O2_inlet_to_Rc_route_still_missing'],
                'Full six-chart operator must keep downstream gap explicit')
            require(old['original_terminal_reference_slope_background_history_P0_join_receipt']==owner.joins['O2_inlet'],
                'Actual reference0-to-O2-slope0 same-function join receipt required')
            require(old['original_terminal_reference_slope_radius_and_phase_seam']=='Rh_reference -> O2_slope',
                'Exact original terminal radius/phase seam required')
            provenance=got['background']['record']['source_provenance']
            require(provenance['chart']=='O2_slope' and ep(provenance['coordinate_box'])==(0,0),
                'Terminal background must be directly queried at the true O2 slope inlet')
            require(got['background']['record']['P0_not_merged_into_pressure_history'],'Terminal O2 P0/P0_Z must stay separate')
            for key in current.RATES:
                same(got['own'][key],got['background']['originals'][key]+got['correction']['values'][key],
                    'Actual O2 slope0 own history mismatch: '+key)
                same(got['own_Z'][key],got['background']['Z_derivatives'][key]+got['correction']['Z_derivatives'][key],
                    'Actual O2 slope0 own Z history mismatch: '+key)
                require(overlaps(got['own'][key],got['cells']['Rh_reference']['own'][key])
                    and overlaps(got['own_Z'][key],got['cells']['Rh_reference']['own_Z'][key]),
                    'Original same-function reference/O2 own cover consistency failed: '+key)
                terminal_rows+=2
            regions[name]=dict(added_complete_original_charts=6,actual_correction_and_own_C0_Z_rows=120,
                actual_incoming_C0_Z_rows=60,serial_composite_cover_overlap_rows=10,route_exit='Rh_reference0 = O2_slope0')
            print('Actual six-chart inlet-to-O2_inlet C1 route checked:',name,flush=True)
        try:owner.route(N=159)
        except ValueError:pass
        else:raise ArithmeticError('Unsafe whole-period frequency accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        actual_added_chart_correction_and_own_C0_Z_rows_checked=rows,
        actual_inherited_C0_Z_rows_checked=incoming_rows,serial_composite_cover_overlap_rows_checked=overlap_rows,
        native_Z_queries_checked=2,complete_added_original_charts_checked=6,regions=regions,
        accepted_original_loop_cutoff_mass_and_nonzero_incoming_affine_theorems_inherited=True,
        original_full_box_signed_sources_and_true_widths_executed=True,
        exact_original_neighbor_radius_and_phase_seams_checked=6,
        accepted_original_background_history_P0_source_joins_bound=7,
        actual_O2_slope0_own_C0_Z_rows_checked=terminal_rows,
        analytic_patch_exp1_background_endpoint_and_exact_unit_log_width_checked=True,
        actual_original_inlet_to_O2_inlet_no_gap_C1_history_covers_checked=True,
        rate0_pressure_C0_Z_memory_and_separate_P0_retained=True,
        bounds_remain_conservative_not_quantitative_terminal_closure=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Middle/O2 inlet C1 PASS:',rows,'actual chart rows;',incoming_rows,'inherited rows',flush=True)
    return result


if __name__=='__main__':run()
