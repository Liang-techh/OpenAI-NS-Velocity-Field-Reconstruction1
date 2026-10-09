"""Whole-Z switch source replay, genuine incoming and own-rate audit."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_switch_finite_N as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source


def encoded(value):return current.current.source.bridge._encode(current.original.serialized(value))


def exact_row(f,actual,record):
    expected=current.decode(f,record)
    assert actual.scale.powers==expected.scale.powers
    assert actual.scale.offset._mpi_==expected.scale.offset._mpi_
    assert actual.coefficient._mpi_==expected.coefficient._mpi_
    assert actual.zero==expected.zero


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():
        assert current.sha(name)==digest,'Source bytes changed: '+name
    counts=dict(whole_Z_cells=0,actual_radial_source_cells=0,fresh_source_queries=0,
        actual_nonlinear_density_C0_Z_rows=0,actual_own_rate_weight_rows=0,
        actual_boundary_C0_Z_rows=0,exact_directed_interval_rows=0,
        independent_original_P0_rows=0,genuine_R100_incoming_rows=0,typed_rejections=0)
    with mp.workdps(540):
        owner=current.WholeZSwitchFiniteN()
        assert owner.identity==saved['source_family'] and owner.N==saved['candidate_N']
        c=owner.c;right=mp.mpf(-1)
        assert saved['original_prefix_replay_binding']==owner.prefix_binding
        for ends,packet in zip(current.current.source.CELLS,saved['source_cells']):
            op=owner.owner(ends);f=op.flow
            assert current.ep(op.background.Z)[0]==right;right=current.ep(op.background.Z)[1]
            assert packet['candidate_N']==owner.N and packet['source_identity']==owner.identity
            assert packet['exact_Z_cell']==list(ends)
            for row,record in zip(op.background.reference.P0,packet['exact_common_P0_axial5']):
                exact_row(f,row,record);counts['independent_original_P0_rows']+=1
            incoming={name:list(rows) for name,rows in op.incoming.items()}
            for name in current.RATES:
                for n in range(2):
                    exact_row(f,incoming[name][n],packet['actual_R100_correction_incoming_C0_Z'][name][n])
                    counts['genuine_R100_incoming_rows']+=1
            assert packet['actual_R100_incoming_not_zeroed_or_replaced_by_background']
            assert packet['higher_background_jets_not_fabricated_from_C0_Z_correction']
            for chart,window in zip(current.CHARTS,packet['actual_source_windows']):
                assert chart==window['chart'] and len(window['source_cells'])==2
                total={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
                for left,right_,cell in zip(current.PARTITION,current.PARTITION[1:],window['source_cells']):
                    query=cell['source'];assert query['candidate_N']==owner.N
                    assert query['exact_Z_cell']==list(ends) and query['chart']==chart
                    assert query['exact_left']==list(left) and query['exact_right']==list(right_)
                    assert query['real_same_N_R100_correction_supplied']
                    assert query['actual_phase_Z_exact_zero'] and query['phase_average_cancellation_not_claimed']
                    for row,record in zip(op.background.reference.P0,query['exact_common_P0_axial5']):
                        exact_row(f,row,record);counts['independent_original_P0_rows']+=1
                    if chart=='first_switch' and left==(0,1):
                        fresh=owner.query(ends,chart,left,right_)
                        same_source(encoded(fresh),query,counts);counts['fresh_source_queries']+=1
                    recovered=current.decode(f,query['original_generic_source'])
                    got=current.decode(f,query['original_primitive_values'])
                    E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
                    density=current.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got,owner.N)
                    same_source(encoded(density),query['original_signed_five_density_C0_Z'],counts)
                    assert recovered['full_inertial_linear_quadratic_pressure_meridional_sectors_retained']
                    logA=query['actual_original_A_log_absolute_upper']
                    if logA is not None:assert current.ep(c.ln(owner.N)-current.read(c,logA))[0]>0
                    for name,rate in current.RATES.items():
                        weights=current.switch.own_weights(op.first,chart,left,right_,rate)
                        mass,decay,tail=weights
                        for row,key in zip(weights,('true_full_mass','true_cell_decay','true_suffix_decay')):
                            exact_row(f,row,cell['own_rate_weights'][name][key])
                        assert cell['own_rate_weights'][name]['own_rate']==str(rate)
                        assert cell['own_rate_weights'][name]['true_radial_measure_applied_once']
                        if name=='p':assert current.ep(decay.coefficient)==(1,1) and current.ep(tail.coefficient)==(1,1)
                        pair=[current.bounds.symmetric(f,current.bounds.magnitude(f,density[part][name])*mass*tail)
                              for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):
                            exact_row(f,row,cell['signed_cell_driver_C0_Z'][name][n]);total[name][n]+=row
                        counts['actual_own_rate_weight_rows']+=3
                    counts['actual_nonlinear_density_C0_Z_rows']+=10;counts['actual_radial_source_cells']+=1
                for name,rate in current.RATES.items():
                    memory=current.switch.own_weights(op.first,chart,(0,1),(1,1),rate)[1]
                    exact_row(f,memory,window['incoming_own_rate_memory'][name])
                    for n in range(2):
                        exact_row(f,incoming[name][n],window['actual_incoming_correction_C0_Z'][name][n])
                        inherited=incoming[name][n]*memory
                        exact_row(f,inherited,window['actual_retained_incoming_C0_Z'][name][n])
                        exact_row(f,total[name][n],window['actual_local_driver_C0_Z'][name][n])
                        incoming[name][n]=inherited+total[name][n]
                        exact_row(f,incoming[name][n],window['actual_exit_correction_C0_Z'][name][n])
                        counts['actual_boundary_C0_Z_rows']+=4
            terminal=owner.query(ends,'post_power',(1,1));counts['fresh_source_queries']+=1
            for name in current.RATES:
                rows=terminal['original_generic_source']['common_own_five_histories_axial5'][name]
                for n in range(2):
                    exact_row(f,incoming[name][n],packet['actual_R110_correction_C0_Z'][name][n])
                    exact_row(f,rows[n],packet['actual_original_R110_background_C0_Z'][name][n])
                    exact_row(f,rows[n]+incoming[name][n],packet['actual_R110_complete_own_history_C0_Z'][name][n])
                    counts['actual_boundary_C0_Z_rows']+=3
            same_source(encoded(terminal['actual_background_source']['fields']),packet['actual_R110_background_field_source'],counts)
            same_source(encoded(terminal['actual_background_source']['histories']),packet['actual_R110_background_six_history_source'],counts)
            same_source(encoded(terminal['original_signed_five_density_C0_Z']['velocities']),packet['actual_R110_candidate_velocity_C0_Z'],counts)
            assert packet['R110_candidate_velocity_is_same_N_full_phase_enclosure']
            exact_row(f,f.scalar(110),packet['exact_R110_radius'])
            assert packet['exact_total_length_identity']=='hb+hb+(log(1.1)-2hb)=log(1.1)'
            # The source-defined interfaces preserve the six background
            # history lists. Equality is inherited from construction, not
            # established by an overlap of two endpoint enclosures.
            assert op.first.inlet_fields is op.background.axis['normalized_actual_fields']
            assert op.first.inlet is op.background.axis['normalized_actual_own_six_moments']
            first_end=op.first.evaluate((1,1))
            assert op.second.inlet is first_end['actual_six_histories']
            assert op.second.inlet_fields is first_end['actual_fields']
            counts['whole_Z_cells']+=1
            print('Whole-Z same-N R110 source/incoming/weight audit: '+str(ends),flush=True)
        assert right==1
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(('-1','-.5'),'bad',(0,1)),
                     lambda:owner.query(('-1','-.5'),'first_switch',(1,1),(0,1))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid source domain must reject')
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        fresh_live_whole_Z_background_functions_and_same_N_density_queries_replayed=True,
        actual_R100_incoming_all_five_C0_Z_drivers_true_weights_and_R110_boundary_checked=True,
        independent_P0_separate_complete_own_history_and_source_defined_joins_checked=True,
        correction_C0_Z_not_padded_into_six_background_jets=True,
        conservative_full_phase_covers_not_phase_average_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole-Z genuine same-N R110 prefix checks passed',flush=True)
    return result


if __name__=='__main__':run()
