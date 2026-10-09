"""Check actual whole-Z slope, scalar reuse, critical q and genuine incoming."""
import copy
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_O2_slope_finite_N as current
import lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N_check as original_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):return current.encode(current.serialized(value))


def independent_transport():
    c=MPIntervalContext();c.dps=180;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    initial={name:[f.scalar('-.2'),f.scalar('.03')] for name in current.RATES}
    local={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
    coefficients=(p.mpf('-.37'),p.mpf('.021'));alpha=p.mpf('.17');weights=0
    for left,right in zip(current.PARTITION,current.PARTITION[1:]):
        l=p.mpf(left[0])/left[1];r=p.mpf(right[0])/right[1]
        for name,rate in current.RATES.items():
            width,suffix,mass,decay,tail=current.slope_weights(f,left,right,rate)
            lam=p.mpf(rate.numerator)/rate.denominator
            for row,value in ((mass,-p.expm1(-lam*(r-l))/lam if lam else r-l),
                    (decay,p.exp(-lam*(r-l))),(tail,p.exp(-lam*(1-r)))):
                lo,hi=current.ep(row.finite_interval(max_log=2000));assert lo<=value<=hi;weights+=1
            for n,value in enumerate(coefficients):
                source=f.scalar(c.mpf(str(value))*c.exp(c.mpf('.17')*c.mpf([l,r])))
                local[name][n]+=current.reference.bounds.symmetric(f,current.reference.bounds.magnitude(f,source)*mass*tail)
    comparisons=0
    for name,rate in current.RATES.items():
        lam=p.mpf(rate.numerator)/rate.denominator
        memory=current.reference.long.kernel_weight(f,c.mpf(1),c.mpf(0),rate)[1]
        for n,value in enumerate(coefficients):
            output=initial[name][n]*memory+local[name][n]
            expected=p.mpf(('-.2','.03')[n])*p.exp(-lam)+value*(p.exp(alpha)-p.exp(-lam))/(lam+alpha)
            lo,hi=current.ep(output.finite_interval(max_log=2000));assert lo<=expected<=hi;comparisons+=1
    return dict(passed=True,independent_nonconstant_negative_C0_nonzero_Z_integrals=comparisons,
        independent_actual_cell_masses_and_decays=weights,actual_nonzero_incoming_memory_retained=True,
        original_total_slope_log_width_one_and_pressure_memory_one=True)


def source_contract(owner,op,source,counts):
    f,c=op.flow,op.c;r=source['original_generic_source'];proof=source['original_full_source_quotients']
    background=source['original_closed_O2_slope_background']
    assert source['exact_common_P0_axial5'] is op.P0 and r['common_original_P0_axial5'] is op.P0
    assert background['raw']['original_P0_normalized_axial5'] is op.P0
    assert background['exact_source_Rm_factor'] is op.Rm_factor
    assert background['original_scalar_partition']==current.SCALAR_CELLS==128
    assert r['original_fixed_source_log_bases'] is f.logs and len(f.logs)==5
    assert r['source_from_same_current_whole_Z_original_O2_slope_background']
    assert source['phase_not_restarted_at_Rref'] and source['actual_phase_Z_exact_zero']
    assert source['phase_average_cancellation_not_claimed'] and source['signed_axis_no_absZ_or_positive_rho_division']
    assert current.ep(source['actual_phase_origin_s_c'])==current.ep(owner.sc)
    assert proof['actual_positive_E']['strict_positive_complete_source_enclosure']
    assert proof['actual_positive_a']['strict_positive_complete_source_enclosure']
    assert all(row.zero for row in proof['actual_b_axial5']) and all(row.zero for row in proof['actual_t0_axial5'])
    assert all(row.zero for row in background['actual_a_axial5'][1:])
    assert source['original_q_C0_Z'][current.Z].zero
    if source['actual_original_A_log_absolute_upper'] is not None:
        assert current.ep(source['actual_original_A_log_absolute_upper'])[1]<current.ep(c.ln(owner.N))[0]
    counts['strict_actual_A_N_budgets']+=1
    for rows in (*r['common_own_five_histories_axial5'].values(),*r['actual_generic_source_numerators'].values(),
            *r['full_signed_inertial_sectors_axial4'].values()):current.reference.parameters.same_source(f,rows)
    for part in ('kernels','Z_derivatives'):
        current.reference.parameters.same_source(f,list(source['original_signed_five_density_C0_Z'][part].values()))
    assert len(r['common_velocity_E_axial5'])==6 and len(r['common_radial_Q_axial4'])==5
    assert not source['actual_high_order_finite_N_correction_jets_installed']


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['full_same_N_correction_prefix_through_O2_slope_exit_installed']
    assert not saved['actual_high_order_finite_N_correction_jets_installed']
    assert all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    symbolic=original_check.symbolic_source();scalar=original_check.independent_scalar_integrals()
    with mp.workdps(240):independent=independent_transport()
    counts=dict(whole_Z_cells=0,closed_slope_source_cells=0,nonlinear_density_C0_Z_rows=0,
        true_own_rate_weight_rows=0,genuine_Rref_incoming_rows=0,slope_exit_correction_rows=0,
        actual_Rref_background_join_rows=0,exact_critical_endpoint_queries=0,strict_actual_A_N_budgets=0,
        exact_directed_interval_rows=0,typed_rejections=0)
    with mp.workdps(540),patch.object(current.original,'OriginalO2SlopeFiniteN',forbidden), \
            patch.object(current.upstream.WholeZRhReferenceFiniteN,'contribution',forbidden), \
            patch.object(current.upstream.WholeZRhReferenceFiniteN,'query',forbidden):
        owner=current.WholeZO2SlopeFiniteN();c=owner.c
        assert saved['source_family']==owner.identity and saved['candidate_N']==owner.N==2**3981
        assert saved['original_source_bindings']==encoded(owner.bindings)
        assert saved['exact_slope_partition']==[list(row) for row in current.PARTITION]
        for ends,stored in zip(current.CELLS,saved['source_cells']):
            live=owner.contribution(ends);same_source(encoded(live),stored,counts)
            op=owner.owner(ends);f=op.flow
            assert live['exact_total_log_length']==1 and live['exact_common_P0_axial5'] is op.P0
            assert op.P0 is op.reference.P0
            join=live['typed_actual_Rref_slope_source_join']
            assert join['exact_original_function_join_by_J_and_masses_zero_at_y0']
            assert join['overlap_only_consistency_not_functional_identity_proof'] and join['same_original_Rref_P0_and_five_history_seed']
            assert join['directed_source_overlap_consistency_rows']==42
            counts['actual_Rref_background_join_rows']+=42;counts['genuine_Rref_incoming_rows']+=10
            assert live['actual_Rref_incoming_binding']['correction_only_not_complete_own_history']
            for cell in live['actual_source_cells']:
                source_contract(owner,op,cell['source'],counts)
                assert current.ep(cell['actual_log_width'])==(mp.mpf('.25'),mp.mpf('.25'))
                pressure=cell['own_rate_weights']['p']
                assert pressure['true_full_mass'].record()==f.scalar('.25').record()
                assert pressure['true_cell_decay'].record()==f.scalar(1).record()
                assert pressure['true_suffix_decay'].record()==f.scalar(1).record()
                counts['closed_slope_source_cells']+=1;counts['nonlinear_density_C0_Z_rows']+=10;counts['true_own_rate_weight_rows']+=15
            for name,rate in current.RATES.items():
                expected=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-c.mpf(rate.numerator)/rate.denominator)
                assert live['incoming_own_rate_memory'][name].record()==expected.record()
                for key in ('actual_Rref_correction_incoming_C0_Z','actual_retained_incoming_C0_Z','actual_local_driver_C0_Z',
                        'actual_O2_slope_exit_correction_C0_Z','actual_original_O2_slope_exit_background_C0_Z',
                        'actual_O2_slope_exit_complete_own_history_C0_Z'):
                    current.reference.parameters.same_source(f,live[key][name])
            end=live['actual_original_O2_slope_exit_source'];source_contract(owner,op,end,counts)
            back=end['original_closed_O2_slope_background']
            assert back['actual_a_axial5'][0].record()==f.scalar(2).record()
            assert current.ep(back['original_J'])==(mp.mpf('.5'),mp.mpf('.5'))
            assert end['original_full_source_quotients']['actual_Delta_axial5'][0].zero
            assert not end['original_q_C0_Z'][current.C0].zero
            assert current.ep(end['original_q_C0_Z'][current.C0].coefficient)[0]>0
            assert end['original_q_C0_Z'][current.Z].zero
            assert end['original_full_source_quotients']['original_q_Z_scope']['original_sigma_identically_one_on_negative_Delta']
            initial=owner.query(ends,(0,1));source_contract(owner,op,initial,counts)
            assert current.ep(initial['original_closed_O2_slope_background']['original_J'])==(0,0)
            assert all(current.ep(row)==(0,0) for row in initial['original_closed_O2_slope_background']['original_dimensionless_masses'])
            counts['exact_critical_endpoint_queries']+=1;counts['slope_exit_correction_rows']+=10;counts['whole_Z_cells']+=1
            print('Whole-Z genuine same-N O2 slope correction audit: '+str(ends),flush=True)
        assert counts['closed_slope_source_cells']==16
        assert owner.scalar_evaluations==saved['original_scalar_cache_evaluations']==6
        assert owner.scalar_hits==saved['original_scalar_cache_hits']==18 and len(owner.scalar_cache)==6
        key=tuple(current.CELLS[0]);good=owner.saved_cells[key]
        for field,value in (('candidate_N',257),('source_identity',{}),('exact_common_P0_axial5',[]),
                ('actual_Rref_correction_C0_Z',{}),('exact_Rref_radius',encoded(owner.owner(key).flow.scalar(1)))):
            bad=copy.deepcopy(good);bad[field]=value;owner.saved_cells[key]=bad;owner.inlets.pop(key,None)
            try:owner.inlet(key)
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Wrong current Rref incoming accepted: '+field)
        owner.saved_cells[key]=good;owner.inlets.pop(key,None)
        foreign=MPIntervalContext();foreign.dps=c.dps
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(key,(-1,1)),lambda:owner.query(key,(2,1)),
                lambda:owner.query(key,(1,1),(0,1)),lambda:owner.scalars(foreign,foreign.mpf('.5'),128)):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid whole-Z slope context/domain accepted')
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        original_whole_Z_slope_history_functions_variable_shear_and_signed_C0_Z_drivers_checked=True,
        genuine_Rref_correction_only_incoming_true_memory_and_slope_exit_checked=True,
        exact_critical_a_two_and_positive_q_squared_eta_over2_retained=True,
        original_scalar_cache_reuses_only_exact_same_context_y_cell_and_partition=True,
        original_scalar_cache_evaluations=owner.scalar_evaluations,original_scalar_cache_hits=owner.scalar_hits,
        independent_slope_source_identities=symbolic,independent_slope_scalar_diagnostics=scalar,
        independent_slope_signed_transport=independent,full_closed_real_axial_and_O2_slope_cover_checked=True,
        actual_high_order_finite_N_correction_jets_installed=False,conservative_full_phase_covers_not_phase_average_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(original_check.__file__).name:current.sha(Path(original_check.__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Whole-Z genuine same-N O2 slope correction checks passed',flush=True);return receipt


if __name__=='__main__':run()
