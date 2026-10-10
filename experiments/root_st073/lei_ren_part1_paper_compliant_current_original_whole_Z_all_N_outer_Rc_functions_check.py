"""Check the changed outer source chain and direct all-N Rc targets."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
from contextlib import ExitStack
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions_check as bridge_check

encoded=lambda value:current.core.current.encode(current.serialized(value))


def forbidden(*args,**kwargs):raise AssertionError('Ancestor integration or fixed-N support/source/target called')


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed prerequisite: '+name)
    with mp.workdps(540):
        owner=current.WholeZAllNOuterRcFunctions()
        assert saved['source_family']==owner.identity and saved['fixed_leading_input_N0']==owner.N0
        assert saved['original_outer_source_binding']==owner.outer_binding
        assert saved['exact_outer_native_partitions']==current.serialized(current.PARTITIONS)
        assert current.PARTITIONS['Rh_reference']==((-5,1),(-4,1),(-3,1),(-2,1),(-1,1),(0,1))
        assert saved['exact_quotient_Z_theorem']==current.core.relative.exact_quotient_theorem()
        counts=dict.fromkeys(current.CHARTS,0);phase_pieces=target_pairs=0
        with ExitStack() as stack:
            for cls in (current.rh.WholeZRhReferenceFiniteN,current.slope.WholeZO2SlopeFiniteN,
                    current.axial.WholeZO2AxialBufferFiniteN,current.o3.WholeZO3RcFiniteN):
                stack.enter_context(patch.object(cls,'query',forbidden))
                stack.enter_context(patch.object(cls,'inlet',forbidden))
            stack.enter_context(patch.object(current.rh.reference.primitives,'all_u_primitive_bounds',forbidden))
            stack.enter_context(patch.object(current.previous.WholeZAllNLongPatchFunctions,'transport',forbidden))
            stack.enter_context(patch.object(current.core.relative,'normalized_targets',forbidden))
            for transport in saved['actual_four_Z_all_N_outer_Rc_transports']:
                ends=tuple(transport['exact_Z_cell']);op=owner.owner.owner(ends);f=op.flow
                decode=lambda row:current.core.relative.restore_row(f,row)
                pairs=lambda rows:{name:current.Pair(*[decode(v) for v in values]) for name,values in rows.items()}
                assert current.previous.equivalent_rows(op.P0,[decode(v) for v in transport['exact_common_P0_axial5']])
                assert current.previous.equivalent_rows([op.Rm_factor,op.Rm_factor*owner.c.exp(1)],
                    [decode(transport['exact_Rm_radius']),decode(transport['exact_Rh_radius'])])
                incoming=owner.normalized_Rh_incoming(ends)
                assert not all(v.value.zero and v.Z.zero for v in incoming.values())
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_Rh_incoming_C0_Z']
                windows=transport['actual_all_N_outer_windows']
                assert [row['actual_chart'] for row in windows]==list(current.CHARTS)
                for window in windows:
                    chart=window['actual_chart'];partition=current.PARTITIONS[chart]
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_incoming_C0_Z']
                    assert [[cell['exact_left'],cell['exact_right']] for cell in window['actual_all_N_source_cells']]==current.serialized(list(zip(partition,partition[1:])))
                    local={name:current.Pair(f.scalar(0),f.scalar(0)) for name in current.RATES}
                    for cell in window['actual_all_N_source_cells']:
                        phase_pieces+=bridge_check.phase_cover(cell['full_phase_inverse_function_alternatives'])
                        assert cell['physical_Jacobian_applied_once'] and cell['all_N_phase_cover_not_fixed_N_phase_reuse'] and cell['bridge_collar_not_inherited']
                        assert current.ep(current.previous.patch_source.read(owner.c,cell['epsilon']))==current.ep(current.core.epsilon_cover(owner.c,owner.N0))
                        left,right=tuple(cell['exact_left']),tuple(cell['exact_right']);density=pairs(cell['N_scaled_density_C0_Z']);contributions={}
                        if chart=='O3_power':
                            assert all(v.value.zero and v.Z.zero for v in density.values())
                            assert all(alt['original_inverse_proof']['geometry']=='flat' for alt in cell['full_phase_inverse_function_alternatives'])
                        for name,rate in current.RATES.items():
                            weight=owner.weights(ends,chart,left,right,rate)
                            assert encoded(weight)==cell['original_own_rate_weights'][name]
                            contribution=current.scale(density[name],weight['original_mass']*weight['suffix_decay'])
                            local[name]=current.add(local[name],contribution);contributions[name]=contribution
                            if name=='p':
                                assert current.ep(weight['cell_decay'].finite_interval())==(1,1)
                                assert current.ep(weight['suffix_decay'].finite_interval())==(1,1)
                        assert encoded({k:current.record(v) for k,v in contributions.items()})==cell['actual_N_scaled_contribution_C0_Z']
                        counts[chart]+=1
                    assert encoded({k:current.record(v) for k,v in local.items()})==window['actual_N_scaled_local_C0_Z']
                    memory={name:owner.weights(ends,chart,partition[0],partition[-1],rate)['cell_decay'] for name,rate in current.RATES.items()}
                    assert encoded(memory)==window['original_incoming_own_rate_memory']
                    assert current.ep(memory['p'].finite_interval())==(1,1)
                    incoming={name:current.add(current.scale(incoming[name],memory[name]),local[name]) for name in current.RATES}
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_outgoing_C0_Z']
                    if chart=='O2_axial':
                        width=owner.weights(ends,chart,partition[0],partition[-1],current.RATES['p'])['physical_log_width']
                        assert width._mpi_==owner.c.expm1(40)._mpi_
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_Rc_correction_C0_Z']
                terminal=owner.leading_packet(ends,'O3_power',(2,1));background=terminal['original_closed_O3_background']
                amplitude=[terminal['original_roots']['E'][part] for part in (current.C0,current.Z)]
                assert encoded(amplitude)==transport['actual_positive_Rc_normalized_amplitude_C0_Z']
                assert current.previous.equivalent_rows([background['original_positive_formal_mu']],
                    [f.factor((0,0,0,0,0),owner.logmu)])
                assert current.core.relative.exact_row(decode(transport['exact_Rc_radius']),background['actual_source_radius'])
                target=current.normalized_targets(f,incoming,amplitude,background['original_positive_formal_mu'],owner.logmu)
                for name,value in target.items():assert encoded(value)==transport[name],name
                target_pairs+=len(target['actual_uniform_N_scaled_target_C0_Z'])
            assert counts==dict(Rh_reference=20,O2_slope=16,O2_axial=16,O2_buffer=12,O3_transition=16,O3_power=8)
            fresh_sources=[];comparisons=0
            for chart in current.CHARTS:
                partition=current.PARTITIONS[chart];left,right=partition[0],partition[1]
                fresh=owner.query(current.CELLS[0],chart,left,right)
                stored=saved['actual_four_Z_all_N_outer_Rc_transports'][0]['actual_all_N_outer_windows'][current.CHARTS.index(chart)]['actual_all_N_source_cells'][0]
                assert encoded({k:current.record(v) for k,v in fresh['N_scaled_density_C0_Z'].items()})==stored['N_scaled_density_C0_Z']
                assert encoded(fresh['source']['full_phase_inverse_function_alternatives'])==stored['full_phase_inverse_function_alternatives']
                packet=fresh['source']['original_source_packet']
                assert packet['fixed_N_density_support_phase_and_inlet_not_called'] and packet['bridge_collar_not_inherited']
                assert packet['exact_common_P0_axial5'] is owner.owner.owner(current.CELLS[0]).P0
                if chart=='O3_power':
                    assert packet['original_q_C0_Z'][current.C0].zero and packet['original_q_C0_Z'][current.Z].zero
                    assert current.ep(packet['original_power_flat_admission'])[0]>=0
                velocity=fresh['source']['original_E_V_C0_Z']
                for alternative in fresh['source']['full_phase_inverse_function_alternatives']:
                    primitive=alternative['primitive_C0_Z_phi']
                    direct=current.core.bridge.SIGNED_DENSITY(velocity['E'].value,velocity['E'].Z,
                        velocity['V'].value,velocity['V'].Z,primitive,owner.N0)
                    normalized=current.core.normalized_drivers(velocity['E'],velocity['V'],primitive,owner.N0,owner.c.mpf(1)/owner.N0)
                    for key,q in normalized['drivers'].items():
                        for part,row in (('kernels',q.value),('Z_derivatives',q.Z)):
                            lo,hi=current.ep((direct[part][key]*owner.N0-row).coefficient)
                            assert lo<=0<=hi,(chart,key,part);comparisons+=1
                fresh_sources.append(chart)
            for chart,left,right in (('Rh_reference',(-6,1),(-5,1)),('O2_axial',(1,1),(0,1)),
                    ('O3_power',(2,1),(3,1)),('actual_patch',(1,1),(2,1))):
                try:owner.leading_packet(current.CELLS[0],chart,left,right)
                except ValueError:pass
                else:raise AssertionError('Invalid outer source domain accepted')
            try:owner.at_N(current.CELLS[0],'O2_axial',(0,1),(1,4),N=owner.N0-1)
            except ValueError:pass
            else:raise AssertionError('Unadmitted frequency accepted')
        for flag in (*current.OPEN,'actual_N_squared_averaging_cancellation_proved','actual_terminal_controls_installed',
                'actual_global_frequency_admitted','physical_original_exterior_five_targets_closed'):
            assert saved[flag] is False
        for flag in ('actual_uniform_N_scaled_Rc_targets_installed','actual_all_regions_all_N_adapter_installed'):
            assert saved[flag] is True
        hashes=dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            actual_uniform_outer_integration_cells_replayed=sum(counts.values()),actual_chart_cell_counts=counts,
            exact_full_phase_pieces_checked=phase_pieces,actual_six_outer_C1_source_cells_recomputed=fresh_sources,
            actual_fixed_N0_same_primitive_density_consistency_comparisons=comparisons,
            actual_uniform_Rc_target_C0_Z_pairs_replayed=target_pairs,
            all_four_Z_live_P0_Rm_Rh_Rc_positive_mu_amplitude_incoming_and_pressure_memory_checked=True,
            original_frontends_pure_parent_recipes_physical_measures_and_quiet_incoming_retained=True,
            unchanged_integrations_fixed_N_outer_queries_inlets_support_and_target_function_forbidden=True,
            actual_uniform_N_scaled_Rc_targets_installed=True,actual_all_regions_all_N_adapter_installed=True,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.core.current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: actual whole-Z full-phase all-N outer functions, Rc transport and direct relative targets',flush=True)
    return receipt


if __name__=='__main__':run()
