"""Check whole-Z repaired reference source and actual same-N Rref incoming."""
import copy
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_reference_finite_N as current
import lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N_check as original_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):
    return current.upstream.upstream.upstream.source.bridge._encode(current.serialized(value))


def independent_transport():
    """Exact nonconstant signed forcing with nonzero correction incoming."""
    c=MPIntervalContext();c.dps=180;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    initial={name:[f.scalar('-.2'),f.scalar('.03')] for name in current.RATES}
    local={name:[f.scalar(0),f.scalar(0)] for name in current.RATES}
    coefficients=(p.mpf('-.37'),p.mpf('.021'));alpha=p.mpf('.17')
    weights=0
    for left,right in zip(current.PARTITION,current.PARTITION[1:]):
        l,r=left[0],right[0]
        for name,rate in current.RATES.items():
            width,suffix,mass,decay,tail=current.reference_weights(f,left,right,rate)
            lam=p.mpf(rate.numerator)/rate.denominator
            for row,value in ((mass,-p.expm1(-lam)/lam if lam else p.mpf(1)),
                    (decay,p.exp(-lam)),(tail,p.exp(lam*r))):
                lo,hi=current.ep(row.finite_interval(max_log=2000));assert lo<=value<=hi;weights+=1
            assert current.ep(width)==(1,1) and current.ep(suffix)==(-r,-r)
            for n,value in enumerate(coefficients):
                source=f.scalar(c.mpf(str(value))*c.exp(c.mpf('.17')*c.mpf([l,r])))
                local[name][n]+=current.reference.bounds.symmetric(f,
                    current.reference.bounds.magnitude(f,source)*mass*tail)
    comparisons=0
    for name,rate in current.RATES.items():
        lam=p.mpf(rate.numerator)/rate.denominator
        memory=current.upstream.long.kernel_weight(f,c.mpf(5),c.mpf(0),rate)[1]
        for n,value in enumerate(coefficients):
            output=initial[name][n]*memory+local[name][n]
            expected=p.mpf(('-.2','.03')[n])*p.exp(-5*lam)+value*(-p.expm1(-5*(lam+alpha)))/(lam+alpha)
            lo,hi=current.ep(output.finite_interval(max_log=2000));assert lo<=expected<=hi;comparisons+=1
    return dict(passed=True,independent_nonconstant_negative_C0_nonzero_Z_integrals=comparisons,
        independent_actual_cell_masses_and_decays=weights,actual_nonzero_incoming_memory_retained=True,
        original_five_unit_log_width_and_pressure_memory_one=True)


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['full_same_N_correction_prefix_through_Rref_installed']
    assert not saved['actual_high_order_finite_N_correction_jets_installed']
    assert all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    symbolic=original_check.symbolic_source()
    with mp.workdps(240):independent=independent_transport()
    counts=dict(whole_Z_cells=0,closed_reference_source_cells=0,nonlinear_density_C0_Z_rows=0,
        true_own_rate_weight_rows=0,genuine_Rh_incoming_rows=0,Rref_correction_rows=0,
        actual_Rh_background_join_rows=0,strict_actual_A_N_budgets=0,
        exact_directed_interval_rows=0,typed_rejections=0)
    with mp.workdps(540),patch.object(current.original,'OriginalRhReferenceFiniteN',forbidden), \
            patch.object(current.upstream.WholeZRmPatchFiniteN,'contribution',forbidden), \
            patch.object(current.upstream.WholeZRmPatchFiniteN,'query',forbidden):
        owner=current.WholeZRhReferenceFiniteN();c=owner.c
        assert saved['source_family']==owner.identity and saved['candidate_N']==owner.N==2**3981
        assert saved['original_source_bindings']==encoded(owner.bindings)
        assert 'exact_reference_phase' not in owner.bindings
        assert saved['exact_reference_partition']==[list(row) for row in current.PARTITION]
        for ends,stored in zip(current.upstream.upstream.upstream.source.CELLS,saved['source_cells']):
            live=owner.contribution(ends);same_source(encoded(live),stored,counts)
            op=owner.owner(ends);f=op.flow
            assert live['exact_total_log_length']==5 and live['exact_common_P0_axial5'] is op.P0
            assert op.P0 is op.reference.P0 and len(f.logs)==5
            join=live['typed_actual_Rh_reference_source_join']
            assert join['source_function_identity_from_original_full_weight_unique_leading_map']
            assert join['overlap_only_consistency_not_functional_identity_proof'] and join['shared_original_P0']
            assert join['directed_source_overlap_consistency_rows']==42
            counts['actual_Rh_background_join_rows']+=42;counts['genuine_Rh_incoming_rows']+=10
            inlet=live['actual_Rh_incoming_binding']
            assert inlet['candidate_N']==owner.N and inlet['correction_only_not_complete_own_history']
            assert inlet['no_old_N257_label_incoming_or_phase_receipt_transplanted']
            for cell in live['actual_source_cells']:
                source=cell['source'];background=source['actual_background_source']
                r=source['original_generic_source'];proof=source['original_full_source_quotients']
                assert source['exact_common_P0_axial5'] is op.P0 and background['original_P0_normalized_axial5'] is op.P0
                assert r['common_original_P0_axial5'] is op.P0 and r['original_fixed_source_log_bases'] is f.logs
                assert r['source_from_same_current_whole_Z_repaired_Rh_background'] and r['exact_E_y_equals_E_over10_and_a_four_fifths']
                assert source['phase_not_restarted_at_Rh'] and source['actual_phase_Z_exact_zero']
                assert source['phase_average_cancellation_not_claimed'] and source['signed_axis_no_absZ_or_positive_rho_division']
                assert current.ep(source['actual_phase_origin_s_c'])==current.ep(owner.sc)
                assert current.ep(cell['actual_log_width'])==(1,1)
                assert proof['source_positive_E']['strict_positive_complete_source_enclosure']
                assert all(row.zero for row in proof['actual_nonzero_b_axial5'])
                assert all(row.zero for row in proof['actual_t0_axial5'])
                assert not source['original_q_C0_Z'][current.C0].zero and source['original_q_C0_Z'][current.Z].zero
                assert current.ep(proof['actual_correlated_a_axial5'][0].coefficient)==current.ep(c.mpf('.8'))
                assert all(row.zero for row in proof['actual_correlated_a_axial5'][1:])
                if source['actual_original_A_log_absolute_upper'] is not None:
                    assert current.ep(source['actual_original_A_log_absolute_upper'])[1]<current.ep(c.ln(owner.N))[0]
                for rows in (*r['common_own_five_histories_axial5'].values(),
                        *r['actual_generic_source_numerators'].values(),*r['full_signed_inertial_sectors_axial4'].values()):
                    current.upstream.parameters.same_source(f,rows)
                assert len(r['common_velocity_E_axial5'])==6 and len(r['common_radial_Q_axial4'])==5
                for part in ('kernels','Z_derivatives'):
                    current.upstream.parameters.same_source(f,list(source['original_signed_five_density_C0_Z'][part].values()))
                pressure=cell['own_rate_weights']['p']
                assert pressure['true_full_mass'].record()==f.scalar(1).record()
                assert pressure['true_cell_decay'].record()==f.scalar(1).record()
                assert pressure['true_suffix_decay'].record()==f.scalar(1).record()
                counts['closed_reference_source_cells']+=1;counts['nonlinear_density_C0_Z_rows']+=10
                counts['true_own_rate_weight_rows']+=15;counts['strict_actual_A_N_budgets']+=1
            for name,rate in current.RATES.items():
                expected=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-5*c.mpf(rate.numerator)/rate.denominator)
                assert live['incoming_own_rate_memory'][name].record()==expected.record()
                for key in ('actual_Rh_correction_incoming_C0_Z','actual_retained_incoming_C0_Z',
                        'actual_local_driver_C0_Z','actual_Rref_correction_C0_Z',
                        'actual_original_Rref_background_C0_Z','actual_Rref_complete_own_history_C0_Z'):
                    current.upstream.parameters.same_source(f,live[key][name])
            assert live['actual_original_Rref_source']['exact_left']==[0,1]
            counts['Rref_correction_rows']+=10;counts['whole_Z_cells']+=1
            print('Whole-Z genuine same-N Rh..Rref correction audit: '+str(ends),flush=True)
        assert counts['closed_reference_source_cells']==20
        key=tuple(current.upstream.upstream.upstream.source.CELLS[0]);good=owner.saved_cells[key]
        for field,value in (('candidate_N',257),('source_identity',{}),('exact_common_P0_axial5',[]),
                ('actual_Rh_correction_C0_Z',{}),('exact_Rh_radius',encoded(owner.owner(key).flow.scalar(1)))):
            bad=copy.deepcopy(good);bad[field]=value;owner.saved_cells[key]=bad;owner.inlets.pop(key,None)
            try:owner.inlet(key)
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Wrong current Rh incoming accepted: '+field)
        owner.saved_cells[key]=good;owner.inlets.pop(key,None)
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(key,(-6,1)),
                lambda:owner.query(key,(0,1),(-1,1)),lambda:owner.query(key,(1,1))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid current whole-Z reference domain accepted')
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        original_whole_Z_repaired_Rh_reference_source_and_signed_C0_Z_drivers_checked=True,
        genuine_Rh_correction_only_incoming_true_memory_and_Rref_exit_checked=True,
        independent_reference_source_identities=symbolic,independent_reference_signed_transport=independent,
        full_closed_real_axial_and_Rh_Rref_reference_cover_checked=True,
        actual_high_order_finite_N_correction_jets_installed=False,
        conservative_full_phase_covers_not_phase_average_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(original_check.__file__).name:current.sha(Path(original_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Whole-Z genuine same-N Rh..Rref correction checks passed',flush=True);return receipt


if __name__=='__main__':run()
