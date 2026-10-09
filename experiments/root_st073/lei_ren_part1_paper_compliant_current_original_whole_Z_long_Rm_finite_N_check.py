"""Audit new whole-Z long/reference functions and genuine same-N Rm transport."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_finite_N as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_whole_Z_switch_finite_N_check import exact_row


def encoded(value):return current.current.source.bridge._encode(current.reference.serialized(value))


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and all(saved[key] is False for key in current.OPEN)
    assert saved['whole_Z_leading_Rm_implicit_controls_solved'] is False
    for name,digest in saved['input_hashes'].items():
        assert current.sha(name)==digest,'Source bytes changed: '+name
    counts=dict(whole_Z_cells=0,actual_radial_source_cells=0,fresh_source_queries=0,
        actual_nonlinear_density_C0_Z_rows=0,actual_own_rate_weight_rows=0,
        actual_boundary_C0_Z_rows=0,exact_directed_interval_rows=0,
        independent_original_P0_rows=0,genuine_R110_incoming_rows=0,
        actual_finite_kernel_backends=0,source_defined_background_joins=0,typed_rejections=0)
    with mp.workdps(540), patch.object(current.endpoint,'OriginalLongReshapeEndpoint',forbidden), \
            patch.object(current.background,'OriginalReferenceRestoreFunctions',forbidden):
        owner=current.WholeZLongRmFiniteN()
        assert owner.identity==saved['source_family'] and owner.N==saved['candidate_N']
        assert saved['original_source_bindings']==encoded(owner.bindings)
        c=owner.c;right=mp.mpf(-1)
        for ends,packet in zip(current.current.source.CELLS,saved['source_cells']):
            op=owner.owner(ends);f=op.flow;ref=op.reference
            assert current.ep(op.Z)[0]==right;right=current.ep(op.Z)[1]
            assert packet['candidate_N']==owner.N and packet['source_identity']==owner.identity
            assert packet['exact_Z_cell']==list(ends)
            for row,record in zip(ref.P0,packet['exact_common_P0_axial5']):
                exact_row(f,row,record);counts['independent_original_P0_rows']+=1
            # Rebuild genuine anchored G derivatives from its rational G' source.
            source=owner.source.background.source
            gradient=source.amplitude.gradient_jets(op.Z,5)
            G=op.original_B_source['actual_anchored_G']
            for n in range(1,6):assert G[n]._mpi_==(gradient[n-1]/n)._mpi_
            same_source(encoded(op.original_B_source),packet['actual_original_anchored_B_source'],counts)
            same_source(encoded(op.terminal),packet['actual_original_finite_long_endpoint'],counts)
            assert ref.initial is op.terminal['actual_terminal_normalized_six_history_shapes']
            assert op.long.histories is op.R110_histories and op.long.phi is op.R110_fields['phi']
            assert ref.P0 is ref.reference((0,1))['original_P0_normalized_axial5']
            assert ref.reference((0,1))['actual_normalized_six_history_shapes'] is ref.initial
            assert ref.restoration((0,1))['actual_centered_histories'] is ref.reference((1,1))['actual_centered_histories']
            assert ref.postrestore((-7,1))['actual_centered_histories'] is ref.restoration((1,1))['actual_centered_histories']
            counts['source_defined_background_joins']+=5
            counts['actual_finite_kernel_backends']+=len(current.kernels.KINDS)
            assert op.terminal['original_nonzero_inlet_memory_and_body_tail_errors_retained']
            incoming={name:list(rows) for name,rows in op.incoming.items()}
            for name in current.RATES:
                for n in range(2):
                    exact_row(f,incoming[name][n],packet['actual_R110_correction_incoming_C0_Z'][name][n])
                    counts['genuine_R110_incoming_rows']+=1
            assert packet['original_full_window_and_nonzero_incoming_memory_retained']
            assert packet['actual_R110_background_not_replaced_by_finite_correction']
            for chart,window in zip(current.CHARTS,packet['actual_source_windows']):
                assert chart==window['chart'] and len(window['source_cells'])==2
                length=owner.T if chart=='long_reshape' else ref.gap if chart=='reference' else c.mpf(1)
                assert length._mpi_==current.read(c,window['actual_window_length'])._mpi_
                total={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
                for left,right_,cell in zip(current.PARTITION,current.PARTITION[1:],window['source_cells']):
                    query=cell['source'];assert query['candidate_N']==owner.N
                    assert query['exact_Z_cell']==list(ends) and query['chart']==chart
                    assert query['exact_left']==list(left) and query['exact_right']==list(right_)
                    assert query['actual_radius_phase_Z_exact_zero'] and query['phase_average_cancellation_not_claimed']
                    assert query['actual_same_N_R110_incoming_available']
                    for row,record in zip(ref.P0,query['exact_common_P0_axial5']):
                        exact_row(f,row,record);counts['independent_original_P0_rows']+=1
                    if chart=='long_reshape' and left==(0,1):
                        fresh=owner.query(ends,chart,left,right_)
                        same_source(encoded(fresh),query,counts);counts['fresh_source_queries']+=1
                    recovered=current.upstream.decode(f,query['original_generic_source'])
                    got=current.upstream.decode(f,query['original_primitive_values'])
                    E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
                    density=current.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got,owner.N)
                    same_source(encoded(density),query['original_signed_five_density_C0_Z'],counts)
                    assert recovered['full_inertial_linear_quadratic_pressure_meridional_sectors_retained']
                    logA=query['actual_original_A_log_absolute_upper']
                    if logA is not None:assert current.ep(c.ln(owner.N)-current.read(c,logA))[0]>0
                    l,r=current.long.fraction(left),current.long.fraction(right_)
                    width=length*c.mpf((r-l).numerator)/(r-l).denominator
                    suffix=length*c.mpf((1-r).numerator)/(1-r).denominator
                    for name,rate in current.RATES.items():
                        mass,decay,tail=current.long.kernel_weight(f,width,suffix,rate)
                        for row,key in zip((mass,decay,tail),('true_full_mass','true_cell_decay','true_suffix_decay')):
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
                    memory=current.long.kernel_weight(f,length,c.mpf(0),rate)[1]
                    exact_row(f,memory,window['incoming_own_rate_memory'][name])
                    for n in range(2):
                        exact_row(f,incoming[name][n],window['actual_incoming_correction_C0_Z'][name][n])
                        inherited=incoming[name][n]*memory
                        exact_row(f,inherited,window['actual_retained_incoming_C0_Z'][name][n])
                        exact_row(f,total[name][n],window['actual_local_driver_C0_Z'][name][n])
                        incoming[name][n]=inherited+total[name][n]
                        exact_row(f,incoming[name][n],window['actual_exit_correction_C0_Z'][name][n])
                        counts['actual_boundary_C0_Z_rows']+=4
                if chart=='long_reshape':
                    for name in current.RATES:
                        for n in range(2):exact_row(f,incoming[name][n],packet['actual_Rsh_correction_C0_Z'][name][n])
            terminal=owner.query(ends,'postrestore',(1,1));counts['fresh_source_queries']+=1
            for name in current.RATES:
                rows=terminal['original_generic_source']['common_own_five_histories_axial5'][name]
                for n in range(2):
                    exact_row(f,incoming[name][n],packet['actual_Rm_correction_C0_Z'][name][n])
                    exact_row(f,rows[n],packet['actual_original_Rm_background_C0_Z'][name][n])
                    exact_row(f,rows[n]+incoming[name][n],packet['actual_Rm_complete_own_history_C0_Z'][name][n])
                    counts['actual_boundary_C0_Z_rows']+=3
            Rm=ref.postrestore((-6,1))
            same_source(encoded(Rm),packet['actual_original_Rm_background_function'],counts)
            assert Rm['original_P0_normalized_axial5'] is ref.P0
            assert Rm['geometry']['exact_reference_offset']==[-6,1]
            exact_row(f,f.factor((0,5,0,0,0),10*owner.logC+c.ln(110)-6),packet['exact_Rm_radius'])
            assert packet['exact_total_length_identity']=='T+gap+1+1=10*(logC+logP)-6'
            assert current.ep(owner.T+ref.gap+2-(10*(owner.logC+ref.logP)-6))[0]<=0<=current.ep(owner.T+ref.gap+2-(10*(owner.logC+ref.logP)-6))[1]
            counts['whole_Z_cells']+=1
            print('Whole-Z same-N actual Rm source/incoming audit: '+str(ends),flush=True)
        assert right==1
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(('-1','-.5'),'bad',(0,1)),
                     lambda:owner.query(('-1','-.5'),'long_reshape',(1,1),(0,1))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid source domain must reject')
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        fresh_live_whole_Z_long_reference_functions_and_same_N_density_queries_replayed=True,
        actual_R110_incoming_all_five_C0_Z_drivers_true_weights_and_Rm_boundary_checked=True,
        independent_P0_separate_complete_own_history_and_source_defined_joins_checked=True,
        genuine_finite_kernels_tail_errors_and_anchored_G_derivative_convention_checked=True,
        correction_C0_Z_not_padded_into_six_background_jets=True,
        conservative_full_phase_covers_not_phase_average_cancellation=True,
        whole_Z_leading_Rm_implicit_controls_solved=False,replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole-Z genuine same-N Rm prefix checks passed',flush=True)
    return result


if __name__=='__main__':run()
