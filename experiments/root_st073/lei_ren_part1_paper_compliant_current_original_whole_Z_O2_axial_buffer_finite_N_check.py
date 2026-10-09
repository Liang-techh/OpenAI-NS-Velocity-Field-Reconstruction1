"""Independent O2 identities, true physical transport and whole-Z Rd replay."""
import copy
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_O2_axial_buffer_finite_N as current
import lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N_check as original_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):return current.encode(current.serialized(value))


def independent_transport():
    c=MPIntervalContext();c.dps=180;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    comparisons=weights=0;alpha=p.mpf('.17');coefficients=(p.mpf('-.37'),p.mpf('.021'))
    # Rescale only this independent small fixture; production remains Md40.
    with patch.object(current.original,'MD',2):
        for chart,partition in current.PARTITIONS.items():
            W=p.expm1(2) if chart=='axial' else p.mpf(11)
            initial={name:[f.scalar('-.2'),f.scalar('.03')] for name in current.RATES}
            local={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
            for left,right in zip(partition,partition[1:]):
                l=p.mpf(left[0])/left[1];r=p.mpf(right[0])/right[1]
                yl=p.exp(2*l)-1 if chart=='axial' else l
                yr=p.exp(2*r)-1 if chart=='axial' else r
                for name,rate in current.RATES.items():
                    width,suffix,mass,decay=current.original.physical_weights(f,chart,left,right,rate)
                    lam=p.mpf(rate.numerator)/rate.denominator
                    for row,value in ((width,yr-yl),(suffix,W-yr)):
                        lo,hi=current.ep(row);assert lo<=value<=hi;weights+=1
                    for row,value in ((mass,-p.expm1(-lam*(yr-yl))/lam if lam else yr-yl),
                            (decay,p.exp(-lam*(W-yr)))):
                        interval=row.finite_interval(max_log=2000) if hasattr(row,'scale') else row
                        lo,hi=current.ep(interval);assert lo<=value<=hi;weights+=1
                    for n,value in enumerate(coefficients):
                        source=f.scalar(c.mpf(str(value))*c.exp(c.mpf('.17')*c.mpf([yl,yr])))
                        local[name][n]+=current.reference.bounds.symmetric(f,current.reference.bounds.magnitude(f,source)*mass*decay)
            output=current.original.formal_affine_transport(f,initial,local,c.mpf(W))
            for name,rate in current.RATES.items():
                lam=p.mpf(rate.numerator)/rate.denominator
                for n,value in enumerate(coefficients):
                    expected=p.mpf(('-.2','.03')[n])*p.exp(-lam*W)+value*(p.exp(alpha*W)-p.exp(-lam*W))/(lam+alpha)
                    lo,hi=current.ep(output[name][n].finite_interval(max_log=2000));assert lo<=expected<=hi;comparisons+=1
    # At the real enormous Md40 length, verify logs and nonzero formal memory.
    actual=0
    for chart,partition in current.PARTITIONS.items():
        W=p.expm1(40) if chart=='axial' else p.mpf(11)
        for left,right in zip(partition,partition[1:]):
            l=p.mpf(left[0])/left[1];r=p.mpf(right[0])/right[1]
            yl=p.exp(40*l) if chart=='axial' else l;yr=p.exp(40*r) if chart=='axial' else r
            end=p.exp(40) if chart=='axial' else p.mpf(11)
            for name,rate in current.RATES.items():
                width,suffix,mass,decay=current.original.physical_weights(f,chart,left,right,rate)
                for row,value in ((width,yr-yl),(suffix,end-yr)):
                    lo,hi=current.ep(row);assert lo<=value<=hi;actual+=1
                current.reference.parameters.positive_source(f,f.scalar(mass),'independent_actual_kernel_mass')
                if rate:
                    log=-p.mpf(rate.numerator)/rate.denominator*(end-yr)
                    lo,hi=current.ep(decay.scale.offset);assert lo<=log<=hi;assert not decay.zero
                else:assert decay.record()==f.scalar(1).record()
        seed={name:[f.scalar(1),f.scalar(1)] for name in current.RATES}
        zero={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
        out=current.original.formal_affine_transport(f,seed,zero,c.mpf(W))
        assert all(not row.zero for pair in out.values() for row in pair)
        assert out['p'][0].record()==f.scalar(1).record()
    return dict(passed=True,independent_nonconstant_negative_C0_nonzero_Z_integrals=comparisons,
        independent_rescaled_mass_width_suffix_decay_checks=weights,independent_actual_Md40_width_suffix_checks=actual,
        actual_huge_length_memory_is_formal_positive_not_underflowed=True,pressure_memory_exactly_one=True)


def source_contract(owner,op,source,counts):
    f,c=op.flow,op.c;r=source['original_generic_source'];proof=source['original_full_source_quotients']
    back=source['original_closed_O2_axial_buffer_background'];a=back['actual_a_axial5']
    assert source['exact_common_P0_axial5'] is r['common_original_P0_axial5'] is op.P0
    assert back['raw']['original_P0_normalized_axial5'] is op.P0 and back['exact_source_Rm_factor'] is op.Rm_factor
    assert r['original_fixed_source_log_bases'] is f.logs and len(f.logs)==5
    assert r['source_from_same_current_whole_Z_original_O2_axial_buffer_background']
    assert source['phase_not_restarted_at_O2_slope_exit'] and source['actual_phase_Z_exact_zero']
    assert source['phase_average_cancellation_not_claimed'] and source['signed_axis_no_absZ_or_positive_rho_division']
    assert current.ep(source['actual_phase_origin_s_c'])==current.ep(owner.sc)
    assert proof['actual_positive_E']['strict_positive_complete_source_enclosure']
    assert proof['actual_positive_a']['strict_positive_complete_source_enclosure']
    assert a[0].record()==f.scalar(2).record() and all(row.zero for row in a[1:])
    assert source['exact_critical_a2_and_collected_Delta_b_squared_over2']
    assert all(not row.zero for row in back['positive_formal_history_decays'].values())
    for rows in (*r['common_own_five_histories_axial5'].values(),*r['actual_generic_source_numerators'].values(),
            *r['full_signed_inertial_sectors_axial4'].values()):current.reference.parameters.same_source(f,rows)
    for part in ('kernels','Z_derivatives'):
        current.reference.parameters.same_source(f,list(source['original_signed_five_density_C0_Z'][part].values()))
    if source['actual_original_A_log_absolute_upper'] is not None:
        assert current.ep(source['actual_original_A_log_absolute_upper'])[1]<current.ep(c.ln(owner.N))[0]
    counts['strict_actual_A_N_budgets']+=1
    if source['actual_chart']=='buffer':
        assert all(current.ep(row)==(0,0) for row in back['original_cutoff_ordinary_y_derivatives'])
        assert all(row.zero for row in r['common_velocity_V_axial5'])
        assert all(row.zero for row in proof['actual_Delta_axial5'])
        assert not source['original_q_C0_Z'][current.C0].zero and source['original_q_C0_Z'][current.Z].zero
        assert current.ep(source['original_q_C0_Z'][current.C0].coefficient)[0]>0
        assert back['original_turnoff_kernels']['buffer_exact_retained_kernel_memory']
        assert current.ep(back['original_turnoff_kernels']['B_mass'])[1]>0
        assert current.ep(back['original_turnoff_kernels']['retained_far_tail'])[1]>0
        counts['positive_critical_buffer_q_with_retained_kernel']+=1
    assert len(r['common_velocity_E_axial5'])==6 and len(r['common_radial_Q_axial4'])==5
    assert not source['actual_high_order_finite_N_correction_jets_installed']


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['full_same_N_correction_prefix_through_Rd_installed']
    assert not saved['actual_high_order_finite_N_correction_jets_installed'] and all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    symbolic=original_check.symbolic_source()
    with mp.workdps(130):scalar=original_check.independent_scalar()
    with mp.workdps(240):independent=independent_transport()
    counts=dict(whole_Z_cells=0,closed_axial_source_cells=0,closed_buffer_source_cells=0,
        nonlinear_density_C0_Z_rows=0,true_own_rate_weight_rows=0,genuine_slope_exit_incoming_rows=0,
        axial_exit_correction_rows=0,Rd_exit_correction_rows=0,functional_join_consistency_rows=0,
        positive_critical_buffer_q_with_retained_kernel=0,strict_actual_A_N_budgets=0,
        exact_directed_interval_rows=0,typed_rejections=0)
    with mp.workdps(540),patch.object(current.original,'OriginalO2AxialBufferFiniteN',forbidden), \
            patch.object(current.upstream.WholeZO2SlopeFiniteN,'contribution',forbidden), \
            patch.object(current.upstream.original,'OriginalO2SlopeFiniteN',forbidden):
        owner=current.WholeZO2AxialBufferFiniteN();c=owner.c
        assert saved['source_family']==owner.identity and saved['candidate_N']==owner.N==2**3981
        assert saved['original_source_bindings']==encoded(owner.bindings)
        assert saved['exact_original_parameter_binding']==encoded(owner.parameters)
        assert saved['exact_chart_partitions']==encoded(current.PARTITIONS)
        for ends,stored in zip(current.CELLS,saved['source_cells']):
            live=owner.contribution(ends);same_source(encoded(live),stored,counts)
            op=owner.owner(ends);f=op.flow
            assert live['exact_common_P0_axial5'] is op.P0 is op.reference.P0
            for join in live['typed_actual_O2_axial_buffer_source_joins'].values():
                assert join['exact_original_function_join_by_closed_kernel_and_decay_identities']
                assert join['overlap_only_consistency_not_functional_identity_proof'] and join['same_original_P0_and_slope_history_seed']
                assert join['directed_source_overlap_consistency_rows']==42
                counts['functional_join_consistency_rows']+=42
            counts['genuine_slope_exit_incoming_rows']+=10
            assert live['actual_slope_exit_incoming_binding']['correction_only_not_complete_own_history']
            assert live['actual_O2_axial_buffer_charts']['buffer']['actual_chart_incoming_C0_Z'] is \
                live['actual_O2_axial_buffer_charts']['axial']['actual_chart_exit_correction_C0_Z']
            for chart,part in live['actual_O2_axial_buffer_charts'].items():
                width=c.expm1(40) if chart=='axial' else c.mpf(11)
                assert current.ep(part['exact_total_log_width']-width)[0]<=0<=current.ep(part['exact_total_log_width']-width)[1]
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
            source_contract(owner,op,owner.query(ends,'axial',(0,1)),counts)
            source_contract(owner,op,owner.query(ends,'buffer',(0,1)),counts)
            assert live['actual_Rd_correction_C0_Z'] is live['actual_O2_axial_buffer_charts']['buffer']['actual_chart_exit_correction_C0_Z']
            assert live['actual_original_Rd_source']['original_closed_O2_axial_buffer_background']['raw']['original_P0_normalized_axial5'] is op.P0
            counts['axial_exit_correction_rows']+=10;counts['Rd_exit_correction_rows']+=10;counts['whole_Z_cells']+=1
            print('Whole-Z genuine same-N O2 axial/buffer correction audit: '+str(ends),flush=True)
        assert counts['closed_axial_source_cells']==16 and counts['closed_buffer_source_cells']==12
        assert owner.kernel_evaluations==saved['original_kernel_cache_evaluations']==6
        assert owner.kernel_hits==saved['original_kernel_cache_hits']==38 and len(owner.kernel_cache)==6
        key=tuple(current.CELLS[0]);good=owner.saved_cells[key]
        for field,value in (('candidate_N',257),('source_identity',{}),('exact_common_P0_axial5',[]),
                ('actual_O2_slope_exit_correction_C0_Z',{}),('exact_O2_slope_exit_radius',encoded(owner.owner(key).flow.scalar(1)))):
            bad=copy.deepcopy(good);bad[field]=value;owner.saved_cells[key]=bad;owner.inlets.pop(key,None)
            try:owner.inlet(key)
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Wrong current slope-exit incoming accepted: '+field)
        owner.saved_cells[key]=good;owner.inlets.pop(key,None)
        foreign=MPIntervalContext();foreign.dps=c.dps
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(key,'axial',(-1,1)),
                lambda:owner.query(key,'axial',(2,1)),lambda:owner.query(key,'axial',(1,1),(0,1)),
                lambda:owner.query(key,'buffer',(12,1)),lambda:owner.query(key,'other',(0,1)),
                lambda:owner.kernels(foreign,foreign.mpf(1),40,128,800)):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid whole-Z O2 context/domain accepted')
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        original_whole_Z_axial_buffer_history_functions_critical_quotients_and_signed_C0_Z_drivers_checked=True,
        genuine_slope_exit_correction_only_incoming_true_physical_memory_and_Rd_exit_checked=True,
        exact_critical_a_two_collected_Delta_and_positive_buffer_q_retained=True,
        original_kernel_cache_reuses_only_exact_same_context_y_Md_partition_window=True,
        original_kernel_cache_evaluations=owner.kernel_evaluations,original_kernel_cache_hits=owner.kernel_hits,
        independent_axial_buffer_source_identities=symbolic,independent_kernel_diagnostics=scalar,
        independent_physical_signed_transport=independent,full_closed_real_axial_and_O2_axial_buffer_cover_checked=True,
        actual_high_order_finite_N_correction_jets_installed=False,conservative_full_phase_covers_not_phase_average_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(original_check.__file__).name:current.sha(Path(original_check.__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Whole-Z genuine same-N O2 axial/buffer correction checks passed',flush=True);return receipt


if __name__=='__main__':run()
