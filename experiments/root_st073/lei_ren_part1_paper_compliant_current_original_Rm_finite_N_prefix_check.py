"""Independent nonzero-inlet transport and typed live current-source binding."""
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rm_finite_N_prefix as current

fields,ep=current.fields,current.ep


def independent_nonzero_transport():
    c=MPIntervalContext();c.dps=150;p=mp.mp.clone();p.dps=100
    macro=current.prefix.macro.macro
    f=macro.fields.MacroFlow(c,c.ln(c.mpf('.002')),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    incoming={};drivers={};truth={};comparisons=0
    for k,rate in current.RATES.items():
        lam=p.mpf(rate.numerator)/rate.denominator
        # Different negative and positive C0/Z data detect reset and mixing.
        initial=(p.mpf('-.37'),p.mpf('.021'));density=(p.mpf('.11'),p.mpf('-.013'))
        incoming[k]=[f.scalar(c.mpf(v)) for v in initial]
        mass=p.mpf(1) if not rate else -p.expm1(-lam)/lam
        drivers[k]=[f.scalar(c.mpf(str(v*mass))) for v in density]
        truth[k]=[initial[n]*p.exp(-lam)+density[n]*mass for n in (0,1)]
    got=current.weighted.terminal.affine_transport(f,incoming,drivers,c.mpf(1))
    for k,pair in got.items():
        for row,value in zip(pair,truth[k]):
            lo,hi=ep(row.finite_interval(max_log=2000));epsilon=p.mpf('1e-90')
            assert lo-epsilon<=value<=hi+epsilon;comparisons+=1
    zero={k:[f.scalar(0),f.scalar(0)] for k in current.RATES}
    quiet=current.weighted.terminal.affine_transport(f,incoming,zero,c.mpf(1))
    for k,pair in quiet.items():
        rate=current.RATES[k];lam=p.mpf(rate.numerator)/rate.denominator
        for n,row in enumerate(pair):
            value=p.mpf(('-.37','.021')[n])*p.exp(-lam)
            lo,hi=ep(row.finite_interval(max_log=2000));epsilon=p.mpf('1e-90')
            assert lo-epsilon<=value<=hi+epsilon and not row.zero;comparisons+=1
    assert current.base.encoded(fields.serialized(quiet['p']))==current.base.encoded(fields.serialized(incoming['p']))
    return dict(passed=True,independent_negative_C0_nonzero_Z_affine_and_quiet_memory_comparisons=comparisons,
        pressure_rate_zero_memory_exactly_one=True,zero_local_source_does_not_reset_incoming=True,
        analytic_reference_not_source_replay=True)


def native(owner,report):
    incoming_rows=driver_rows=outgoing_rows=cells=phase_points=0
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(fields.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        assert packet['exact_source_Rm_factor'] is op.Rm_factor
        assert packet['full_original_logarithmic_interval_width']._mpi_==c.mpf(1)._mpi_
        saved=owner.saved[current.prefix.NAME]['frames'][label]
        local=owner.saved[current.weighted.NAME]['frames'][label]
        baseline=owner.saved[current.weighted.previous.previous.NAME]['frames'][label]
        counts=[0,0,0]
        for i,key in enumerate(('actual_current_Rm_incoming_correction_C0_Z',
            'accepted_complete_local_Rm_Rh_driver_C0_Z','actual_current_Rh_correction_C0_Z')):
            assert set(packet[key])==set(current.RATES)
            for pair in packet[key].values():
                for row in pair:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;counts[i]+=1
        incoming_rows+=counts[0];driver_rows+=counts[1];outgoing_rows+=counts[2]
        assert current.base.encoded(fields.serialized(packet['actual_current_Rm_incoming_correction_C0_Z']))==saved['actual_current_Rm_incoming_correction_C0_Z']
        for k in current.RATES:
            actual=packet['accepted_complete_local_Rm_Rh_driver_C0_Z'][k]
            assert current.base.encoded(fields.serialized(actual))==local['actual_averaged_whole_patch_local_integral_C0_Z'][k]
            assert current.base.encoded(fields.serialized(actual[1]))==baseline['actual_whole_patch_local_defect_integral_C0_Z'][k][1]
            assert not any(row.zero for row in packet['actual_current_Rm_incoming_correction_C0_Z'][k])
        # Independently form e^-lambda * actual incoming + accepted local.
        # This checks C0/Z ordering and own-rate assignment on real rows.
        for k,rate in current.RATES.items():
            lam=c.mpf(rate.numerator)/rate.denominator
            memory=c.mpf(1) if not rate else c.exp(-lam)
            assert packet['original_incoming_Rm_to_Rh_decays'][k]._mpi_==memory._mpi_
            for n in (0,1):
                with mp.workdps(c.dps+40):
                    rebuilt=packet['actual_current_Rm_incoming_correction_C0_Z'][k][n]*memory+packet['accepted_complete_local_Rm_Rh_driver_C0_Z'][k][n]
                assert current.base.encoded(fields.serialized(rebuilt))==current.base.encoded(fields.serialized(packet['actual_current_Rh_correction_C0_Z'][k][n]))
        assert ep(packet['original_incoming_Rm_to_Rh_decays']['p'])==(1,1)
        assert packet['bound_complete_local_source_cell_count']==10
        assert packet['real_finite_N_Rm_incoming_correction_is_supplied']
        assert packet['genuine_current_finite_N_prefix_through_Rh_boundary_enclosures_supplied']
        assert packet['no_ancestor_producer_checker_or_source_integration_rerun']
        assert not packet['ordinary_Z_Nminus2_averaging_installed']
        assert all(packet[k] is False for k in fields.previous.OPEN)
        cells+=packet['bound_complete_local_source_cell_count'];phase_points+=len(current.weighted.PARTITION)
        print('Real source-bound incoming transported Rm-to-Rh checked',label,flush=True)
    label='0';op=owner.owner(label)
    source=owner.saved[current.prefix.NAME]['frames'][label]
    local=owner.saved[current.weighted.NAME]['frames'][label]
    mapper=current.weighted.previous.phase.RmRadiusPhase(op,owner.family,owner.upstream.parameters.upstream.upstream.upstream.repair)
    phases={x:current.base.encoded(fields.serialized(mapper.point(x,257))) for x in current.weighted.PARTITION}
    rejects=0
    for key,value in (('source_family','bad'),('source_frame','.5'),('candidate_N',1024),('exact_common_P0_axial5',[]),
        ('genuine_current_finite_N_R0_R100_R110_Rm_boundary_enclosures_supplied',False),
        ('local_drivers_composed_with_real_current_prefix_not_arbitrary_zero',False),
        ('actual_current_Rm_incoming_correction_C0_Z',{})):
        bad=copy.deepcopy(source);bad[key]=value
        try:current.guard_prefix(owner.family,label,257,op.P0,bad)
        except ValueError:rejects+=1
        else:raise AssertionError('Invalid actual prefix accepted: '+key)
    mutations=[lambda b:b.__setitem__('source_frame','.5'),lambda b:b.__setitem__('candidate_N',1024),
        lambda b:b.__setitem__('exact_partition',[[1,1],'Rh']),
        lambda b:b.__setitem__('actual_source_weighted_cells',b['actual_source_weighted_cells'][:-1]),
        lambda b:b.__setitem__('actual_averaged_whole_patch_local_integral_C0_Z',{}),
        lambda b:b['actual_source_weighted_cells'][0]['source']['exact_source_Rm_factor'].__setitem__('exact_zero',True),
        lambda b:b['actual_source_weighted_cells'][0]['source'].__setitem__('source_geometry',{}),
        lambda b:b['actual_source_weighted_cells'][0]['source']['actual_endpoint_phases']['left'].__setitem__('exact_original_radius_phase','point midpoint'),
        lambda b:b['actual_source_weighted_cells'][0]['source'].__setitem__('original_full_density_mean_not_zeroed',False),
        lambda b:b.__setitem__('ordinary_Z_Nminus2_averaging_installed',True)]
    for mutation in mutations:
        bad=copy.deepcopy(local);mutation(bad)
        try:current.guard_local(owner.family,label,257,op.P0,op.Rm_factor,bad,phases)
        except ValueError:rejects+=1
        else:raise AssertionError('Invalid local source packet accepted')
    for label,N in (('0',1024),('axis',257)):
        try:owner.contribution(label,N)
        except ValueError:rejects+=1
        else:raise AssertionError('Unsupplied source frame/frequency accepted')
    return dict(passed=True,actual_Rm_incoming_C0_Z_rows=incoming_rows,accepted_complete_local_C0_Z_rows=driver_rows,
        actual_Rh_outgoing_C0_Z_rows=outgoing_rows,complete_local_weighted_source_cells=cells,
        exact_current_endpoint_phase_bindings=phase_points,typed_invalid_source_and_packet_rejections=rejects,
        real_Rm_incoming_and_ordinary_Z_local_source_bitwise_preserved=True,
        source_bound_two_frames_not_whole_Z_functional_closure_or_global_N=True)


def run():
    began=time.monotonic()
    with mp.workdps(190):independent=independent_nonzero_transport()
    print('Independent nonzero source and quiet own-rate memory PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmFiniteNPrefix(require_checked=False)
    actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_nonzero_incoming_reference=independent,actual_live_source=actual,input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.base.encoded(result),indent=2)+'\n').encode())
    print('Current original genuine source-bound finite-N prefix through Rh PASS',flush=True);return result


if __name__=='__main__':run()
