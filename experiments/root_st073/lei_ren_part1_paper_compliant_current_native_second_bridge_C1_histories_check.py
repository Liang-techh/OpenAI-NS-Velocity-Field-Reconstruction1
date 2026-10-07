"""Focused actual phase1-memory, whole second-bridge and phase2 C1 checks.

The original scalar-loop proof and nonlinear-coordinate affine references
are inherited from accepted receipts; this checker targets the new route.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_second_bridge_C1_histories as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;packets=current.packets
require=checks.require;ep=current.ep


def same(left,right,message):
    require(packets.encode(left.record())==packets.encode(right.record()),message)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2,'Actual whole second-bridge C1 records required')
    require(not any(saved.get(k) for k in packets.OPEN),'Second bridge cannot complete global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed second-bridge prerequisite: '+name)
    first_saved=json.loads((HERE/current.prior_stage.NAME).read_bytes())
    inherited=json.loads((HERE/current.prior_stage.RECEIPT).read_bytes())
    require(inherited['all_passed'] and inherited[current.prior_stage.GATE],'Original first-bridge theorem/reference receipt required')
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        p=current.prior_stage
        c1=current.density.NativeDensityC1LocalIntegrals(p.first.NativePhaseFirstJets(p.slow.NativeQSlowJets(
            p.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))))
        owner=current.NativeSecondBridgeC1Histories(p.NativeFirstBridgeC1Histories(
            current.transfer.NativeTrueChartC1Transfer(current.history.NativeC1HistoryTransfer(c1))))
        regions={};rows=0;inherited_rows=0
        for name,old in saved['actual_whole_second_bridge_C1_history_records'].items():
            got=owner.second_bridge(packets.interval(owner.ctx,old['Z_box']),saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual second-bridge C1 record changed: '+name)
            require(old['known_original_first_bridge_phase1_history']==first_saved['actual_whole_first_bridge_C1_history_records'][name],
                'Actual incoming must be the accepted full first-bridge result')
            require(old['no_gap_between_original_inlet_and_phase2'] and old['correction_incoming_is_actual_phase1_not_own_history_or_zero_reset'],
                'Second bridge must inherit actual inlet history with no reset')
            geometry=got['geometry'];require(ep(geometry['coordinate'])==(1,2),'Whole native second-bridge domain required')
            require(not geometry['width'].zero and ep(geometry['width'].coefficient)[0]>0,'Original microscopic hb width must remain positive')
            require(geometry['record']['width_and_endpoints_independent_of_Z'],'Fixed-Z differentiation needs original Z-independent geometry')
            require(geometry['record']['exact_original_radius_Jacobian_identities']['passed'],'Original adjacent radius identities required')
            require(old['original_whole_domain_a_positive_certificate']['source_function_positivity_not_inferred_from_saved_box'],
                'Native source sample cannot supply analytic a lower')
            require(old['source_C0_Z_queried_on_full_native_coordinate_box'],'Whole cell signed-source functions required')
            incoming=got['incoming']['correction']
            require(any(not v.zero for v in incoming['values'].values()),'Actual nonzero first-bridge correction memory was erased')
            require(any(not v.zero for v in incoming['Z_derivatives'].values()),'Actual first-Z correction memory was erased')
            for key in current.RATES:
                decay=got['factors'][key]['decay']
                same(got['correction']['values'][key],decay*incoming['values'][key]+got['contributions'][key],
                    'Actual inherited C0 memory mismatch: '+key)
                same(got['correction']['Z_derivatives'][key],decay*incoming['Z_derivatives'][key]+got['Z_derivatives'][key],
                    'Actual inherited Z memory mismatch: '+key)
                same(got['own'][key],got['background']['originals'][key]+got['correction']['values'][key],
                    'Original phase2 background double-added or omitted: '+key)
                same(got['own_Z'][key],got['background']['Z_derivatives'][key]+got['correction']['Z_derivatives'][key],
                    'Original phase2 background Z memory mismatch: '+key)
            require(ep(got['factors']['p']['decay'].coefficient)==(1,1) and ep(got['factors']['p']['decay'].scale.evaluate())==(0,0),
                'Rate-zero pressure memory must have exact unit transfer coefficient')
            require(old['original_phase2_background_and_separate_P0_Z']['P0_not_merged_into_pressure_history'],
                'Original analytic P0/P0_Z datum must stay separate')
            for value in [*got['own'].values(),*got['own_Z'].values(),*got['correction']['values'].values(),*got['correction']['Z_derivatives'].values()]:
                require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,
                    'Actual phase2 C0/Z rows must share the original common basis/ledger')
            require(not old['global_inlet_to_Rc_histories_admitted'] and not any(old.get(k) for k in packets.OPEN),
                'Conservative two-bridge covers cannot imply terminal/global admission')
            rows+=20;inherited_rows+=10
            regions[name]=dict(actual_phase2_correction_and_own_C0_Z_rows=20,inherited_actual_phase1_C0_Z_rows=10,
                full_second_bridge_coordinate_box=[1,2],whole_route_scope='original sc/2 -> phase1 -> phase2')
            print('Actual full second-bridge phase1 memory and phase2 C1 checked:',name,flush=True)
        try:owner.second_bridge(N=159)
        except ValueError:pass
        else:raise ArithmeticError('Unsafe whole-period candidate frequency accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        actual_phase2_correction_and_own_C0_Z_rows_checked=rows,actual_nonzero_inherited_phase1_C0_Z_rows_checked=inherited_rows,
        native_Z_queries_checked=2,regions=regions,
        inherited_original_scalar_loop_and_cutoff_C1_reference_receipt=current.prior_stage.RECEIPT,
        inherited_scalar_loop_comparisons=inherited['independent_original_loop_and_cutoff_C1_checks']['independent_original_loop_whole_period_value_and_implicit_Z_comparisons'],
        inherited_scalar_cutoff_comparisons=inherited['independent_original_loop_and_cutoff_C1_checks']['independent_original_cutoff_q_C0_Z_union_comparisons'],
        original_true_radius_geometry_and_nonzero_incoming_affine_theorems_inherited=True,
        actual_original_inlet_to_phase2_no_gap_C1_history_covers_checked=True,
        rate0_pressure_C0_Z_memory_and_separate_P0_retained=True,
        covers_are_conservative_and_do_not_prove_tight_error_or_terminal_closure=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Whole second-bridge C1 PASS:',rows,'actual phase2 rows;',inherited_rows,'actual inherited rows',flush=True)
    return result


if __name__=='__main__':run()
