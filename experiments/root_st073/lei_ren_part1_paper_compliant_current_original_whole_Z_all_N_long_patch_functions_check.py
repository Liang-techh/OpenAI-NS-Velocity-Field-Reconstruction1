"""Check the changed uniform long/patch chain without ancestor integrations."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_long_patch_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions_check as bridge_check

encoded=lambda value:current.core.current.encode(current.serialized(value))


def forbidden(*args,**kwargs):raise AssertionError('Unchanged integration or fixed-N source/inlet called')


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed prerequisite: '+name)
    with mp.workdps(540):
        owner=current.WholeZAllNLongPatchFunctions()
        assert saved['source_family']==owner.identity and saved['fixed_leading_input_N0']==owner.N0
        assert saved['original_long_patch_source_binding']==owner.long_patch_binding
        assert saved['exact_long_native_partition']==[[0,1],[1,2],[1,1]]
        assert saved['exact_patch_partition']==current.serialized(current.PATCH_PARTITION)
        cells=long_cells=patch_cells=phase_pieces=0
        with patch.object(current.previous.WholeZAllNSwitchFunctions,'transport',forbidden), \
                patch.object(current.core.WholeZAllNBridgeSourceFunctions,'transport',forbidden), \
                patch.object(current.long_signed.source_module.WholeZLongRmFiniteN,'query',forbidden), \
                patch.object(current.patch_source.WholeZRmPatchFiniteN,'query',forbidden), \
                patch.object(current.patch_source.WholeZRmPatchFiniteN,'inlet',forbidden), \
                patch.object(current.core.current.WholeZSharpBridgeFunctions,'query',forbidden), \
                patch.object(current.core.current.WholeZSharpBridgeFunctions,'local_integral',forbidden):
            for transport in saved['actual_four_Z_all_N_long_patch_transports']:
                ends=tuple(transport['exact_Z_cell']);op=owner.owner.owner(ends);f=op.flow
                decode=lambda row:current.core.relative.restore_row(f,row)
                pairs=lambda rows:{name:current.Pair(*[decode(v) for v in values]) for name,values in rows.items()}
                assert current.equivalent_rows(op.P0,[decode(row) for row in transport['exact_common_P0_axial5']])
                assert current.equivalent_rows([op.Rm_factor],[decode(transport['exact_Rm_radius'])])
                assert current.equivalent_rows([f.scalar(110)],
                    [decode(owner.R110_rows[ends]['exact_R110_radius'])])
                incoming=owner.normalized_R110_incoming(ends)
                assert not all(q.value.zero and q.Z.zero for q in incoming.values())
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_R110_incoming_C0_Z']
                windows=transport['actual_all_N_long_patch_windows']
                assert [row['actual_chart'] for row in windows]==list(current.CHARTS)
                for window in windows:
                    chart=window['actual_chart'];length=owner.window_length(ends,chart)
                    assert encoded(length)==window['actual_original_window_length']
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_incoming_C0_Z']
                    partition=current.PATCH_PARTITION if chart=='actual_patch' else current.LONG_PARTITION
                    assert [[cell['exact_left'],cell['exact_right']] for cell in window['actual_all_N_source_cells']]==current.serialized(list(zip(partition,partition[1:])))
                    local={name:current.Pair(f.scalar(0),f.scalar(0)) for name in current.RATES}
                    for cell in window['actual_all_N_source_cells']:
                        phase_pieces+=bridge_check.phase_cover(cell['full_phase_inverse_function_alternatives'])
                        assert cell['physical_Jacobian_applied_once'] and cell['all_N_phase_cover_not_fixed_N_phase_reuse'] and cell['bridge_collar_not_inherited']
                        assert current.ep(current.patch_source.read(owner.c,cell['epsilon']))==current.ep(current.core.epsilon_cover(owner.c,owner.N0))
                        left=tuple(cell['exact_left']);right='Rh' if cell['exact_right']=='Rh' else tuple(cell['exact_right'])
                        density=pairs(cell['N_scaled_density_C0_Z']);contributions={}
                        for name,rate in current.RATES.items():
                            weight=owner.weights(ends,chart,left,right,rate)
                            assert current.ep(weight['physical_log_width'])[0]>0 and current.ep(weight['original_mass'].coefficient)[0]>0
                            contribution=current.scale(density[name],weight['original_mass']*weight['suffix_decay'])
                            local[name]=current.add(local[name],contribution);contributions[name]=contribution
                            assert encoded(weight)==cell['original_own_rate_weights'][name]
                            if name=='p':
                                assert current.ep(weight['cell_decay'].finite_interval())==(1,1)
                                assert current.ep(weight['suffix_decay'].finite_interval())==(1,1)
                        assert encoded({k:current.record(v) for k,v in contributions.items()})==cell['actual_N_scaled_contribution_C0_Z']
                        cells+=1
                        if chart=='actual_patch':patch_cells+=1
                        else:long_cells+=1
                    assert encoded({k:current.record(v) for k,v in local.items()})==window['actual_N_scaled_local_C0_Z']
                    memory={name:current.patch_source.long.kernel_weight(f,length,owner.c.mpf(0),rate)[1] for name,rate in current.RATES.items()}
                    assert encoded(memory)==window['original_incoming_own_rate_memory']
                    assert current.ep(memory['p'].finite_interval())==(1,1)
                    incoming={name:current.add(current.scale(incoming[name],memory[name]),local[name]) for name in current.RATES}
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_outgoing_C0_Z']
                    if chart=='postrestore':
                        assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_Rm_correction_C0_Z']
                        reference=owner.owner.long_owner(ends).reference
                        terminal=current.long_signed.source_module.reference.background_cell(reference,'postrestore',owner.c.mpf(1))
                        assert current.equivalent_rows([op.Rm_factor],[terminal['actual_physical_radius']])
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_Rh_correction_C0_Z']
                assert transport['normalized_correction_only_incoming_background_P0_distinct'] and transport['fixed_N_local_densities_inlets_and_phase_not_used']
            assert (long_cells,patch_cells,cells)==(32,44,76)
            fresh_sources=[];comparisons=0
            for chart,left,right,window_index,cell_index in (('long_reshape',(0,1),(1,2),0,0),
                    ('actual_patch',(49,40),(5,4),4,1)):
                fresh=owner.query(current.CELLS[0],chart,left,right)
                stored=saved['actual_four_Z_all_N_long_patch_transports'][0]['actual_all_N_long_patch_windows'][window_index]['actual_all_N_source_cells'][cell_index]
                assert encoded({k:current.record(v) for k,v in fresh['N_scaled_density_C0_Z'].items()})==stored['N_scaled_density_C0_Z']
                assert encoded(fresh['source']['full_phase_inverse_function_alternatives'])==stored['full_phase_inverse_function_alternatives']
                packet=fresh['source']['original_source_packet']
                assert packet['fixed_N_density_and_phase_and_inlet_not_called'] and packet['bridge_collar_not_inherited']
                assert packet['exact_common_P0_axial5'] is owner.owner.owner(current.CELLS[0]).P0
                velocity=fresh['source']['original_E_V_C0_Z']
                for alternative in fresh['source']['full_phase_inverse_function_alternatives']:
                    primitive=alternative['primitive_C0_Z_phi']
                    direct=current.core.bridge.SIGNED_DENSITY(velocity['E'].value,velocity['E'].Z,
                        velocity['V'].value,velocity['V'].Z,primitive,owner.N0)
                    normalized=current.core.normalized_drivers(velocity['E'],velocity['V'],primitive,owner.N0,owner.c.mpf(1)/owner.N0)
                    for key,q in normalized['drivers'].items():
                        for part,row in (('kernels',q.value),('Z_derivatives',q.Z)):
                            lo,hi=current.ep((direct[part][key]*owner.N0-row).coefficient)
                            assert lo<=0<=hi,(chart,key,part)
                            comparisons+=1
                fresh_sources.append(chart)
            for chart,left,right in (('first_switch',(0,1),(1,2)),('long_reshape',(1,1),(0,1)),('actual_patch',(1,2),(1,1))):
                try:owner.leading_packet(current.CELLS[0],chart,left,right)
                except ValueError:pass
                else:raise AssertionError('Invalid long/patch source domain accepted')
            try:owner.at_N(current.CELLS[0],'long_reshape',(0,1),(1,2),N=owner.N0-1)
            except ValueError:pass
            else:raise AssertionError('Unadmitted frequency accepted')
        for flag in (*current.OPEN,'actual_uniform_N_scaled_Rc_targets_installed','actual_all_regions_all_N_adapter_installed',
                'actual_N_squared_averaging_cancellation_proved','actual_terminal_controls_installed',
                'actual_global_frequency_admitted','physical_original_exterior_five_targets_closed'):
            assert saved[flag] is False
        hashes=dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            actual_uniform_long_integration_cells_replayed=long_cells,actual_uniform_patch_integration_cells_replayed=patch_cells,
            actual_uniform_long_patch_integration_cells_replayed=cells,exact_full_phase_pieces_checked=phase_pieces,
            actual_long_and_patch_C1_source_cells_recomputed=fresh_sources,
            actual_fixed_N0_same_primitive_density_consistency_comparisons=comparisons,
            all_four_Z_live_P0_Rm_radius_R110_incoming_positive_widths_pressure_memory_retained=True,
            exact_source_T_gap_original_long_partition_and_symbolic_Rh_retained=True,
            unchanged_ancestor_integrations_fixed_N_queries_and_patch_inlet_forbidden=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,actual_all_regions_all_N_adapter_installed=False,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.core.current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: actual whole-Z full-phase all-N long/patch functions and normalized Rh transport',flush=True)
    return receipt


if __name__=='__main__':run()
