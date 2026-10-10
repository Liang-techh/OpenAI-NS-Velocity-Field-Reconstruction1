"""Replay changed all-N switch measures, nonzero inlet and C1 memories."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_switch_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions_check as bridge_check

encoded=lambda value:current.previous.current.encode(current.serialized(value))


def forbidden(*args,**kwargs):raise AssertionError('Fixed-N source query or unchanged bridge integration called')


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed prerequisite: '+name)
    with mp.workdps(540):
        owner=current.WholeZAllNSwitchFunctions()
        assert saved['source_family']==owner.identity and saved['fixed_leading_input_N0']==owner.N0
        assert saved['original_switch_source_binding']==owner.switch_binding
        assert saved['original_prefix_replay_binding']==owner.prefix_binding
        cells=phase_pieces=0
        with patch.object(current.previous.WholeZAllNBridgeSourceFunctions,'transport',forbidden), \
                patch.object(current.signed_switch.switch.WholeZSwitchFiniteN,'query',forbidden), \
                patch.object(current.previous.current.WholeZSharpBridgeFunctions,'query',forbidden), \
                patch.object(current.previous.current.WholeZSharpBridgeFunctions,'local_integral',forbidden):
            for transport in saved['actual_four_Z_all_N_switch_transports']:
                ends=tuple(transport['exact_Z_cell']);op=owner.owner.owner(ends);native=owner.owner.switch_owner(ends);f=op.flow
                decode=lambda value:current.previous.relative.restore_row(f,value)
                pairs=lambda rows:{name:current.Pair(*[decode(v) for v in pair]) for name,pair in rows.items()}
                assert len(op.P0)==len(transport['exact_common_P0_axial5'])
                assert all(current.previous.relative.exact_row(a,decode(b)) for a,b in zip(op.P0,transport['exact_common_P0_axial5']))
                incoming=owner.normalized_incoming(ends)
                assert not all(row.value.zero and row.Z.zero for row in incoming.values())
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_R100_incoming_C0_Z']
                windows=transport['actual_all_N_switch_windows']
                assert [window['actual_chart'] for window in windows]==list(current.CHARTS)
                for window in windows:
                    chart=window['actual_chart']
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_incoming_C0_Z']
                    local={name:current.Pair(f.scalar(0),f.scalar(0)) for name in current.RATES}
                    assert [(tuple(cell['exact_left']),tuple(cell['exact_right'])) for cell in window['actual_all_N_source_cells']]==list(zip(current.PARTITION,current.PARTITION[1:]))
                    for cell in window['actual_all_N_source_cells']:
                        phase_pieces+=bridge_check.phase_cover(cell['full_phase_inverse_function_alternatives'])
                        assert cell['physical_Jacobian_applied_once'] and cell['all_N_phase_cover_not_fixed_N_phase_reuse'] and cell['first_micro_collar_not_inherited']
                        assert current.ep(current.signed_switch.switch.read(owner.c,cell['epsilon']))==current.ep(current.previous.epsilon_cover(owner.c,owner.N0))
                        length=f.h if chart!='post_power' else f.scalar(owner.c.ln(owner.c.mpf(11)/10))-f.h*2
                        assert current.previous.relative.exact_row(length,decode(cell['original_positive_chart_log_length']))
                        assert current.ep(length.record()['log_absolute_upper'])[1]>-mp.inf
                        density=pairs(cell['N_scaled_density_C0_Z']);contributions={}
                        for name,rate in current.RATES.items():
                            mass,decay,suffix=current.original.own_weights(native.first,chart,tuple(cell['exact_left']),tuple(cell['exact_right']),rate)
                            assert current.ep(mass.coefficient)[0]>0
                            contribution=current.scale(density[name],mass*suffix)
                            contributions[name]=contribution;local[name]=current.add(local[name],contribution)
                            assert encoded(dict(original_mass=mass,cell_decay=decay,suffix_decay=suffix))==cell['original_own_rate_weights'][name]
                        assert encoded({k:current.record(v) for k,v in contributions.items()})==cell['actual_N_scaled_contribution_C0_Z']
                        cells+=1
                    assert encoded({k:current.record(v) for k,v in local.items()})==window['actual_N_scaled_local_C0_Z']
                    memory={name:current.original.own_weights(native.first,chart,(0,1),(1,1),rate)[1] for name,rate in current.RATES.items()}
                    assert encoded(memory)==window['original_incoming_own_rate_memory']
                    assert current.ep(memory['p'].finite_interval())==(1,1)
                    incoming={name:current.add(current.scale(incoming[name],memory[name]),local[name]) for name in current.RATES}
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_outgoing_C0_Z']
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_R110_correction_C0_Z']
                assert transport['actual_nonzero_normalized_R100_incoming_not_reset'] and transport['normalized_incoming_not_passed_to_raw_correction_inlet']
            assert cells==24
            fresh=owner.query(current.CELLS[0],'first_switch',(0,1),(1,2))
            stored=saved['actual_four_Z_all_N_switch_transports'][0]['actual_all_N_switch_windows'][0]['actual_all_N_source_cells'][0]
            assert encoded({k:current.record(v) for k,v in fresh['N_scaled_density_C0_Z'].items()})==stored['N_scaled_density_C0_Z']
            assert encoded(fresh['source']['full_phase_inverse_function_alternatives'])==stored['full_phase_inverse_function_alternatives']
            packet=fresh['source']['original_source_packet']
            assert packet['fixed_N_density_and_fixed_N_phase_provider_not_called'] and packet['first_micro_collar_not_inherited']
            assert packet['exact_common_P0_axial5'] is owner.owner.owner(current.CELLS[0]).P0
            velocity=fresh['source']['original_E_V_C0_Z'];comparisons=0
            for alternative in fresh['source']['full_phase_inverse_function_alternatives']:
                primitive=alternative['primitive_C0_Z_phi']
                direct=current.previous.bridge.SIGNED_DENSITY(velocity['E'].value,velocity['E'].Z,
                    velocity['V'].value,velocity['V'].Z,primitive,owner.N0)
                normalized=current.previous.normalized_drivers(velocity['E'],velocity['V'],primitive,owner.N0,owner.c.mpf(1)/owner.N0)
                for key,q in normalized['drivers'].items():
                    for part,row in (('kernels',q.value),('Z_derivatives',q.Z)):
                        difference=direct[part][key]*owner.N0-row
                        lo,hi=current.ep(difference.coefficient)
                        assert lo<=0<=hi,(key,part)
                        comparisons+=1
            for chart,left,right in (('first_micro',(0,1),(1,2)),('post_power',(-1,1),(1,2)),('first_switch',(1,1),(0,1))):
                try:owner.leading_packet(current.CELLS[0],chart,left,right)
                except ValueError:pass
                else:raise AssertionError('Invalid switch chart or coordinate accepted')
            try:owner.at_N(current.CELLS[0],'first_switch',(0,1),(1,2),N=owner.N0-1)
            except ValueError:pass
            else:raise AssertionError('Unadmitted frequency accepted')
            key=tuple(current.CELLS[0]);original_row=owner.incoming_rows[key]
            zero=owner.owner.owner(key).flow.scalar(0)
            invalid_rows=(dict(original_row,source_identity={}),
                dict(original_row,exact_common_P0_axial5=[]),
                dict(original_row,actual_uniform_N_scaled_R100_correction_C0_Z=
                    encoded({name:current.record(current.Pair(zero,zero)) for name in current.RATES})))
            try:
                for invalid in invalid_rows:
                    owner.incoming_rows[key]=invalid
                    try:owner.normalized_incoming(key)
                    except ValueError:pass
                    else:raise AssertionError('Changed family/P0 or zero-reset incoming accepted')
            finally:owner.incoming_rows[key]=original_row
        for flag in (*current.OPEN,'actual_uniform_N_scaled_Rc_targets_installed',
                'actual_all_regions_all_N_adapter_installed','actual_N_squared_averaging_cancellation_proved',
                'actual_terminal_controls_installed','actual_global_frequency_admitted','physical_original_exterior_five_targets_closed'):
            assert saved[flag] is False
        hashes=dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            actual_uniform_switch_integration_cells_replayed=cells,exact_full_phase_pieces_checked=phase_pieces,
            actual_first_switch_uniform_C1_source_cell_recomputed=True,
            actual_fixed_N0_same_primitive_density_consistency_comparisons=comparisons,
            all_four_Z_original_P0_nonzero_normalized_incoming_positive_widths_pressure_memory_retained=True,
            original_switch_source_AST_binding_checked=True,unchanged_bridge_integrations_and_fixed_N_queries_forbidden=True,
            changed_family_missing_P0_zero_reset_incoming_rejected=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,actual_all_regions_all_N_adapter_installed=False,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.previous.current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: actual whole-Z full-phase all-N switch functions and normalized R110 transport',flush=True)
    return receipt


if __name__=='__main__':run()
