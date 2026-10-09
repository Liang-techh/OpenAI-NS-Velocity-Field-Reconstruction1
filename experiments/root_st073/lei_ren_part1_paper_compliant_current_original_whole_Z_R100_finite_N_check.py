"""Whole-prefix common N, true weights, five density and boundary audit."""
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_R100_finite_N as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source

ep=current.ep


def decode(f,row):
    if isinstance(row,dict):
        if 'formal_positive_scale' in row:return current.switch.first.endpoint.restore_row(f,row)
        if 'lower_exact_mpf_tuple' in row:return current.read(f.c,row)
        return {key:decode(f,value) for key,value in row.items()}
    if isinstance(row,list):return [decode(f,value) for value in row]
    return row


def exact_row(f,actual,record):
    expected=decode(f,record)
    assert actual.scale.powers==expected.scale.powers
    assert actual.scale.offset._mpi_==expected.scale.offset._mpi_
    assert actual.coefficient._mpi_==expected.coefficient._mpi_
    assert actual.zero==expected.zero


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and all(saved[key] is False for key in current.OPEN)
    N=current.phase.candidate_N(saved['candidate_N'])
    selection=saved['current_whole_prefix_N_selection']
    assert N==2**selection['exact_dyadic_power'] and selection['candidate_N']==N
    assert selection['global_all_region_N_cone_or_terminal_repair_admission'] is False
    assert len(selection['whole_prefix_A_requirements'])==40
    for name,digest in saved['input_hashes'].items():
        assert current.sha(name)==digest,'Source input changed: '+name
    counts=dict(whole_Z_cells=0, actual_radial_source_cells=0, actual_nonlinear_density_C0_Z_rows=0,
        actual_own_rate_weight_rows=0, actual_boundary_C0_Z_rows=0, exact_directed_interval_rows=0,
        original_source_owned_zero_inlet_rows=0, common_exponent_requirements=0,
        outward_A_cover_source_rows=0, independent_original_P0_rows=0)
    with mp.workdps(540):
        owner=current.WholeZR100FiniteN()
        assert owner.identity==saved['source_family']
        c=owner.c
        logAmax=-mp.inf
        for record in selection['whole_prefix_A_requirements']:
            logA=record['log_A_absolute_upper']
            if logA is not None:
                upper=ep(current.read(c,logA))[1]
                assert ep(c.ln(N)-c.mpf(upper))[0]>0
                logAmax=max(logAmax,upper)
            counts['common_exponent_requirements']+=1
        if logAmax!=-mp.inf:
            assert ep(current.read(c,selection['maximum_log_A_absolute_upper']))==(logAmax,logAmax)
            assert N==2**max(9,int(mp.ceil(logAmax/mp.log(2)))+1)
        requirements=iter(selection['whole_prefix_A_requirements'])
        right=mp.mpf(-1)
        for ends,packet in zip(current.current.source.CELLS,saved['source_cells']):
            # Fresh leading core/bridge functions establish the exact formal
            # bases. Weight/density audits do not rerun the ancestor chain.
            f,proof=owner.background.source.owner(ends)
            expected_P0=f.jet(current.current.source.IntervalTaylor(c,
                proof['original_pressure_Z0_through_Z6'][:6]))
            for row,record in zip(expected_P0,packet['exact_common_P0_axial5']):
                exact_row(f,row,record);counts['independent_original_P0_rows']+=1
            assert ep(proof['Z'])[0]==right
            right=ep(proof['Z'])[1]
            series=current.current.first.FirstSwitchFunctions.__new__(current.current.first.FirstSwitchFunctions)
            series.flow,series.c=f,c
            assert packet['candidate_N']==N and packet['source_identity']==owner.identity
            assert packet['exact_Z_cell']==list(ends)
            assert all(row['exact_zero'] for rows in packet['actual_current_inlet_correction_C0_Z'].values() for row in rows)
            assert packet['actual_source_owned_initial_condition']['exact_zero_correction_follows_from_defining_Duhamel_lower_bound']
            incoming={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
            counts['original_source_owned_zero_inlet_rows']+=10
            for window,(chart,partition,end) in zip(packet['actual_source_windows'],owner.partitions()):
                assert window['chart']==chart
                total={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
                assert len(window['source_cells'])==len(partition)-1
                for cell,left,right_ in zip(window['source_cells'],partition,partition[1:]):
                    query=cell['source'];assert query['candidate_N']==N
                    for row,record in zip(expected_P0,query['common_original_P0']):
                        exact_row(f,row,record);counts['independent_original_P0_rows']+=1
                    assert query['actual_phase_definition']=='frac(N*log(R/r_minus))'
                    assert query['actual_phase_Z_exact_zero'] and query['phase_average_cancellation_not_claimed']
                    r=decode(f,query['full_signed_generic_source'])
                    got=decode(f,query['original_primitive_values'])
                    requirement=next(requirements)
                    assert requirement['exact_Z_cell']==list(ends) and requirement['chart']==chart
                    primitive_proof=query['original_primitive_proof']
                    if primitive_proof.get('bounded_A_C0_cover_only_not_source_or_A_Z_replacement'):
                        original_A=decode(f,primitive_proof['original_A_C0_before_outward_exponent_cover'])
                        fresh=current.bounded_exponent_cover(f,
                            dict(values={**got,'A':original_A},record={}),N)
                        exact_row(f,fresh['values']['A'],query['original_primitive_values']['A'])
                        assert fresh['values']['A_Z'] is got['A_Z']
                        assert primitive_proof['materialized_A_cap_below_same_explicit_candidate_N']
                        counts['outward_A_cover_source_rows']+=1
                    else:original_A=got['A']
                    assert original_A.zero==requirement['exact_A_zero']
                    if not original_A.zero:
                        upper=ep(original_A.record()['log_absolute_upper'])[1]
                        assert ep(current.read(c,requirement['log_A_absolute_upper']))==(upper,upper)
                    E,V=r['common_velocity_E_axial5'],r['common_velocity_V_axial5']
                    density=current.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got,N)
                    same_source(current.current.source.bridge._encode(current.original.serialized(density)),
                                query['original_signed_five_density_C0_Z'],counts)
                    assert r['full_inertial_linear_quadratic_pressure_meridional_sectors_retained']
                    assert len(r['common_radial_Q_axial4'])==5
                    for name,rate in current.RATES.items():
                        weights=(current.macro.weights(series,left,right_,rate) if chart=='frozen_macro' else
                                 current.original.weights(series,left,right_,end,rate))
                        mass,decay,tail=weights
                        record=cell['own_rate_weights'][name]
                        for row,key in zip(weights,('true_full_mass','true_cell_decay','true_suffix_decay')):
                            exact_row(f,row,record[key])
                            assert row.scale.bases is f.logs and row.ledger is f.ledger
                        assert record['own_rate']==str(rate) and record['true_radial_measure_applied_once']
                        if name=='p':assert ep(decay.coefficient)==(1,1) and ep(tail.coefficient)==(1,1)
                        pair=[current.bounds.symmetric(f,current.bounds.magnitude(f,density[part][name])*mass*tail)
                              for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):
                            exact_row(f,row,cell['signed_cell_driver_C0_Z'][name][n]);total[name][n]+=row
                        counts['actual_own_rate_weight_rows']+=3
                    counts['actual_nonlinear_density_C0_Z_rows']+=10
                    counts['actual_radial_source_cells']+=1
                for name,rate in current.RATES.items():
                    memory=(current.macro.weights(series,partition[0],partition[-1],rate)[1]
                            if chart=='frozen_macro' else
                            current.original.weights(series,partition[0],partition[-1],end,rate)[1])
                    exact_row(f,memory,window['incoming_own_rate_memory'][name])
                    for n in range(2):
                        exact_row(f,incoming[name][n],window['actual_incoming_correction_C0_Z'][name][n])
                        inherited=incoming[name][n]*memory
                        exact_row(f,inherited,window['actual_retained_incoming_C0_Z'][name][n])
                        exact_row(f,total[name][n],window['actual_local_driver_C0_Z'][name][n])
                        incoming[name][n]=inherited+total[name][n]
                        exact_row(f,incoming[name][n],window['actual_exit_correction_C0_Z'][name][n])
                        counts['actual_boundary_C0_Z_rows']+=4
            for name in current.RATES:
                background=decode(f,packet['actual_original_R100_background_C0_Z'][name])
                for n in range(2):
                    exact_row(f,incoming[name][n],packet['actual_current_R100_correction_C0_Z'][name][n])
                    exact_row(f,background[n]+incoming[name][n],packet['actual_current_R100_complete_own_history_C0_Z'][name][n])
                    counts['actual_boundary_C0_Z_rows']+=2
            radius=decode(f,packet['exact_R100_radius'])
            assert ep(radius.coefficient)==(100,100) and radius.scale.powers==(0,0,0,0,0)
            assert packet['actual_micro_macro_source_join']['actual_micro_histories_are_the_live_macro_incoming']
            assert packet['candidate_N_is_not_global_N_or_cone_admission']
            counts['whole_Z_cells']+=1
            print('Whole-Z finite-N density/weight/boundary audit: '+str(ends),flush=True)
        assert right==1
        assert next(requirements,None) is None
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=N,
        exact_whole_Z_source_partition_and_current_common_N_checked=True,
        actual_all_five_density_C0_Z_true_weights_and_affine_boundary_rows_replayed=True,
        complete_own_background_plus_correction_and_separate_P0_checked=True,
        conservative_full_period_enclosures_not_averaged_phase_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
                     Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole-Z genuine current finite-N prefix checks passed',flush=True)
    return result


if __name__=='__main__':run()
