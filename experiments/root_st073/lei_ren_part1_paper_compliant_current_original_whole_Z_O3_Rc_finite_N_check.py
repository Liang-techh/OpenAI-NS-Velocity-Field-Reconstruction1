"""Independent O3 source/mu identities, true transport and whole-Z Rc replay."""
import copy
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N as current
import lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N_check as original_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):return current.encode(current.serialized(value))


def independent_transport():
    c=MPIntervalContext();c.dps=180;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    initial={name:[f.scalar('-.2'),f.scalar('.03')] for name in current.RATES}
    coefficients=(p.mpf('-.37'),p.mpf('.021'));alpha=p.mpf('.17');weights=comparisons=0
    for chart,partition in current.PARTITIONS.items():
        W=p.mpf(current.LENGTHS[chart]);local={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
        for left,right in zip(partition,partition[1:]):
            l=p.mpf(left[0])/left[1];r=p.mpf(right[0])/right[1]
            for name,rate in current.RATES.items():
                width,suffix,mass,decay=current.physical_weights(f,chart,left,right,rate)
                lam=p.mpf(rate.numerator)/rate.denominator
                for row,value in ((width,r-l),(suffix,W-r),(mass,-p.expm1(-lam*(r-l))/lam if lam else r-l),
                        (decay,p.exp(-lam*(W-r)))):
                    interval=row.finite_interval(max_log=2000) if hasattr(row,'scale') else row
                    lo,hi=current.ep(interval);assert lo<=value<=hi;weights+=1
                for n,value in enumerate(coefficients):
                    source=f.scalar(c.mpf(str(value))*c.exp(c.mpf('.17')*c.mpf([l,r])))
                    local[name][n]+=current.reference.bounds.symmetric(f,current.reference.bounds.magnitude(f,source)*mass*decay)
        out=current.original.previous.formal_affine_transport(f,initial,local,c.mpf(W))
        for name,rate in current.RATES.items():
            lam=p.mpf(rate.numerator)/rate.denominator
            for n,value in enumerate(coefficients):
                expected=p.mpf(('-.2','.03')[n])*p.exp(-lam*W)+value*(p.exp(alpha*W)-p.exp(-lam*W))/(lam+alpha)
                lo,hi=current.ep(out[name][n].finite_interval(max_log=2000));assert lo<=expected<=hi;comparisons+=1
    zero={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
    quiet=current.original.previous.formal_affine_transport(f,initial,zero,c.mpf(2))
    assert quiet['p'][0].record()==initial['p'][0].record() and quiet['p'][1].record()==initial['p'][1].record()
    assert all(not row.zero for pair in quiet.values() for row in pair)
    return dict(passed=True,independent_nonconstant_negative_C0_nonzero_Z_integrals=comparisons,
        independent_true_mass_width_suffix_decay_checks=weights,
        quiet_power_keeps_nonzero_real_incoming_with_pressure_memory_one=True)


def source_contract(owner,op,source,counts):
    f,c=op.flow,op.c;r=source['original_generic_source'];proof=source['original_full_source_quotients']
    back=source['original_closed_O3_background'];mu=back['original_positive_formal_mu']
    assert source['exact_common_P0_axial5'] is r['common_original_P0_axial5'] is op.P0
    assert back['raw']['original_P0_normalized_axial5'] is op.P0 and back['exact_source_Rm_factor'] is op.Rm_factor
    assert r['original_fixed_source_log_bases'] is f.logs and len(f.logs)==5
    assert r['source_from_same_current_whole_Z_original_O3_background']
    assert source['phase_not_restarted_at_Rd_or_Rw'] and source['actual_phase_Z_exact_zero']
    assert source['phase_average_cancellation_not_claimed'] and source['signed_axis_no_absZ_or_positive_rho_division']
    assert current.ep(source['actual_phase_origin_s_c'])==current.ep(owner.sc)
    assert proof['actual_positive_E']['strict_positive_complete_source_enclosure']
    assert proof['actual_positive_a']['strict_positive_complete_source_enclosure']
    assert not mu.zero and current.ep(mu.coefficient)==(1,1)
    assert current.ep(mu.scale.offset-owner.logmu)[0]<=0<=current.ep(mu.scale.offset-owner.logmu)[1]
    assert current.ep(back['actual_mu_arithmetic_enclosure_only'])[0]==0
    assert current.ep(back['actual_mu_arithmetic_enclosure_only'])[1]>0
    assert back['exact_source_parameter_not_selected_from_cover']
    Delta=back['actual_collected_Delta_axial5']
    assert current.ep(Delta[0].coefficient)[0]>=0 and all(row.zero for row in Delta[1:])
    assert all(row.zero for row in proof['actual_b_axial5']) and all(row.zero for row in proof['actual_t0_axial5'])
    assert all(row.zero for row in r['common_velocity_V_axial5']) and source['original_q_C0_Z'][current.Z].zero
    for rows in (*r['common_own_five_histories_axial5'].values(),*r['actual_generic_source_numerators'].values(),
            *r['full_signed_inertial_sectors_axial4'].values()):current.reference.parameters.same_source(f,rows)
    for part in ('kernels','Z_derivatives'):
        current.reference.parameters.same_source(f,list(source['original_signed_five_density_C0_Z'][part].values()))
    if source['actual_original_A_log_absolute_upper'] is not None:
        assert current.ep(source['actual_original_A_log_absolute_upper'])[1]<current.ep(c.ln(owner.N))[0]
    counts['strict_actual_A_N_budgets']+=1
    branch=proof['original_q_C0']['branch']
    if branch=='flat':
        assert all(row.zero for row in source['original_q_C0_Z'].values())
        assert all(row.zero for row in source['original_primitive_values'].values())
        assert all(row.zero for part in ('kernels','Z_derivatives') for row in source['original_signed_five_density_C0_Z'][part].values())
        counts['proved_flat_source_queries']+=1
    elif not Delta[0].zero:
        cap=proof['original_q_C0']['original_nonnegative_Delta_q_cap_intersection']
        assert cap['active_source_q_squared_upper']=='eta/2' and cap['flat_source_exact_q_zero']
        assert cap['positive_mu_not_replaced_by_arithmetic_cover']
        counts['source_conditioned_active_flat_union_queries']+=1
    else:
        assert not source['original_q_C0_Z'][current.C0].zero
        assert current.ep(source['original_q_C0_Z'][current.C0].coefficient)[0]>0
        counts['positive_critical_transition_inlet_queries']+=1
    if source['actual_chart']=='power':
        assert branch=='flat' and current.ep(source['actual_same_source_power_flat_log_margin'])[0]>0
        assert not Delta[0].zero and back['original_transition_or_power_kernels']['original_positive_exprel_theta_identity']
    assert not source['actual_high_order_finite_N_correction_jets_installed']


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['full_same_N_correction_prefix_through_Rc_installed']
    assert not saved['actual_high_order_finite_N_correction_jets_installed'] and all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    symbolic=original_check.symbolic_source()
    with mp.workdps(150):scalar=original_check.independent_scalar()
    with mp.workdps(240):independent=independent_transport()
    counts=dict(whole_Z_cells=0,closed_transition_source_cells=0,closed_power_source_cells=0,
        nonlinear_density_C0_Z_rows=0,true_own_rate_weight_rows=0,genuine_Rd_incoming_rows=0,
        Rw_exit_correction_rows=0,Rc_exit_correction_rows=0,functional_join_consistency_rows=0,
        proved_flat_source_queries=0,source_conditioned_active_flat_union_queries=0,
        positive_critical_transition_inlet_queries=0,strict_actual_A_N_budgets=0,
        exact_directed_interval_rows=0,typed_rejections=0)
    with mp.workdps(540),patch.object(current.original,'OriginalO3RcFiniteN',forbidden), \
            patch.object(current.upstream.WholeZO2AxialBufferFiniteN,'contribution',forbidden), \
            patch.object(current.upstream.original,'OriginalO2AxialBufferFiniteN',forbidden):
        owner=current.WholeZO3RcFiniteN();c=owner.c
        assert saved['source_family']==owner.identity and saved['candidate_N']==owner.N==2**3981
        assert saved['original_source_bindings']==encoded(owner.bindings)
        assert saved['exact_original_parameter_binding']==encoded(owner.parameters)
        assert saved['exact_chart_partitions']==encoded(current.PARTITIONS)
        expected_log=c.ln(c.mpf('.001'))-4*(c.exp(40)+11)
        assert current.ep(owner.logmu-expected_log)[0]<=0<=current.ep(owner.logmu-expected_log)[1]
        assert current.ep(owner.logmu)[1]<-1000 and current.ep(owner.power_margin)[0]>0
        assert current.ep(owner.Tw+60*expected_log)[0]<=0<=current.ep(owner.Tw+60*expected_log)[1]
        for ends,stored in zip(current.CELLS,saved['source_cells']):
            live=owner.contribution(ends);same_source(encoded(live),stored,counts)
            op=owner.owner(ends);f=op.flow
            assert live['exact_common_P0_axial5'] is op.P0 is op.reference.P0
            for join in live['typed_actual_O3_source_joins'].values():
                assert join['exact_function_join_by_zero_local_integral_and_decay_one']
                assert join['overlap_only_consistency_not_functional_identity_proof'] and join['same_original_P0_and_five_history_seed']
                assert join['directed_source_overlap_consistency_rows']==42
                counts['functional_join_consistency_rows']+=42
            counts['genuine_Rd_incoming_rows']+=10
            assert live['actual_Rd_incoming_binding']['correction_only_not_complete_own_history']
            transition,power=live['actual_O3_charts']['transition'],live['actual_O3_charts']['power']
            assert power['actual_chart_incoming_C0_Z'] is transition['actual_chart_exit_correction_C0_Z']
            assert all(row.zero for pair in power['actual_chart_local_driver_C0_Z'].values() for row in pair)
            assert any(not row.zero for pair in power['actual_chart_exit_correction_C0_Z'].values() for row in pair)
            for chart,part in live['actual_O3_charts'].items():
                width=c.mpf(current.LENGTHS[chart]);assert current.ep(part['exact_total_log_width'])==current.ep(width)
                for cell in part['actual_source_cells']:
                    source=cell['source'];source_contract(owner,op,source,counts)
                    pressure=cell['own_rate_weights']['p']
                    assert current.ep(pressure['true_full_mass'])==current.ep(cell['actual_log_width'])
                    assert pressure['true_suffix_decay'].record()==f.scalar(1).record()
                    counts['closed_'+chart+'_source_cells']+=1;counts['nonlinear_density_C0_Z_rows']+=10;counts['true_own_rate_weight_rows']+=10
                for name,rate in current.RATES.items():
                    expected=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-width*c.mpf(rate.numerator)/rate.denominator)
                    assert part['incoming_own_rate_memory'][name].record()==expected.record()
                    for key in ('actual_chart_incoming_C0_Z','actual_retained_incoming_C0_Z','actual_chart_local_driver_C0_Z',
                            'actual_chart_exit_correction_C0_Z','actual_original_chart_exit_background_C0_Z',
                            'actual_chart_exit_complete_own_history_C0_Z'):
                        current.reference.parameters.same_source(f,part[key][name])
                source_contract(owner,op,part['actual_original_chart_exit_source'],counts)
            for inlet_chart in ('transition','power'):source_contract(owner,op,owner.query(ends,inlet_chart,(0,1)),counts)
            assert live['actual_Rc_correction_C0_Z'] is power['actual_chart_exit_correction_C0_Z']
            assert live['actual_original_Rc_source']['original_closed_O3_background']['raw']['original_P0_normalized_axial5'] is op.P0
            assert all(current.upstream.exact_row_equal(a,b) for a,b in zip(power['actual_chart_incoming_C0_Z']['p'],power['actual_chart_exit_correction_C0_Z']['p']))
            counts['Rw_exit_correction_rows']+=10;counts['Rc_exit_correction_rows']+=10;counts['whole_Z_cells']+=1
            print('Whole-Z genuine same-N O3 Rc correction audit: '+str(ends),flush=True)
        assert counts['closed_transition_source_cells']==16 and counts['closed_power_source_cells']==8
        assert counts['positive_critical_transition_inlet_queries']==4
        assert owner.kernel_evaluations==saved['original_kernel_cache_evaluations']==6
        assert owner.kernel_hits==saved['original_kernel_cache_hits']==18 and len(owner.kernel_cache)==6
        key=tuple(current.CELLS[0]);good=owner.saved_cells[key]
        for field,value in (('candidate_N',257),('source_identity',{}),('exact_common_P0_axial5',[]),
                ('actual_Rd_correction_C0_Z',{}),('exact_Rd_radius',encoded(owner.owner(key).flow.scalar(1)))):
            bad=copy.deepcopy(good);bad[field]=value;owner.saved_cells[key]=bad;owner.inlets.pop(key,None)
            try:owner.inlet(key)
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Wrong current Rd incoming accepted: '+field)
        owner.saved_cells[key]=good;owner.inlets.pop(key,None)
        foreign=MPIntervalContext();foreign.dps=c.dps;f=owner.owner(key).flow
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(key,'transition',(-1,1)),
                lambda:owner.query(key,'transition',(2,1)),lambda:owner.query(key,'transition',(1,1),(0,1)),
                lambda:owner.query(key,'power',(3,1)),lambda:owner.query(key,'other',(0,1)),
                lambda:owner.kernels(foreign,foreign.mpf(0),foreign.mpf(0)),
                lambda:current.original.nonnegative_Delta_q(f.scalar(2),f.scalar(-1),owner.eta_log,c.ln(2)),
                lambda:current.source_mu_admission(c,c.mpf(0))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid whole-Z O3 context/domain or quiet admission accepted')
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        original_whole_Z_O3_history_functions_positive_mu_critical_quotients_and_signed_C0_Z_drivers_checked=True,
        genuine_Rd_correction_only_incoming_true_memory_and_Rc_exit_checked=True,
        positive_mu_source_Delta_active_transition_and_exact_quiet_power_retained=True,
        original_kernel_cache_reuses_only_exact_same_context_t_mu_partition=True,
        original_kernel_cache_evaluations=owner.kernel_evaluations,original_kernel_cache_hits=owner.kernel_hits,
        independent_O3_source_identities=symbolic,independent_kernel_diagnostics=scalar,
        independent_physical_signed_transport=independent,full_closed_real_axial_and_O3_Rc_cover_checked=True,
        actual_high_order_finite_N_correction_jets_installed=False,conservative_full_phase_covers_not_phase_average_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(original_check.__file__).name:current.sha(Path(original_check.__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Whole-Z genuine same-N O3 Rc correction checks passed',flush=True);return receipt


if __name__=='__main__':run()
