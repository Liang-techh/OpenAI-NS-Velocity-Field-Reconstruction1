"""Contract the actual Rc F pair, then enclose its remaining lower endpoint.

The immutable aggregate endpoint caps are retained as provenance. Their
positive bounds are never subtracted. Other volume and viscous terms keep
the exact accepted quiet179 representation.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import load,sha,digest,require
from .original_quiet_scalar_m15 import materialize_quiet179_view,INDICES,digest as quiet_digest

OPERATION='same_source_Rc_F_pair_replaced_by_remaining_lower_endpoint_enclosure'
ACTIVE_INDICES=[i for i in INDICES if i!=179]


def _interval(iv,tuples):
    return iv.mpf([mp.make_mpf(tuple(v)) for v in tuples])


def _unit_row():
    return dict(approximation='0.0',absolute_error_upper='1',
        raw_directed_error_exact_mpf_tuple=[0,1,0,1],
        source_enclosure_exact_mpf_tuples=[[1,1,0,1],[0,1,0,1]],
        error_applies_to_serialized_decimal=True,source_value_selected=False)


def lower_endpoint_term(sample,source_digest):
    """A symmetric *bound*, with a common directed upper log for C0/C1."""
    iv=MPIntervalContext();iv.dps=250
    values=[_interval(iv,row['weighted_log_bound_exact_mpf_tuples']) for row in sample['components']]
    top=max((v._mpi_[1] for v in values),key=lambda t:mp.make_mpf(t))
    # This atom is a positive error scale only. Its exact upper tuple is used;
    # no exponential, source value, or midpoint is evaluated.
    ctx=mp.mp.clone();ctx.dps=250
    display=ctx.nstr(ctx.make_mpf(top),80)
    display_error=abs(iv.mpf(display)-iv.mpf(mp.make_mpf(top)));error_top=display_error._mpi_[1]
    error_text='0' if error_top[1]==0 else ctx.nstr(ctx.make_mpf(error_top)*ctx.mpf('1.000001'),45)
    logatom=dict(approximation=display,absolute_error_upper=error_text,
        raw_directed_error_exact_mpf_tuple=list(error_top),error_applies_to_serialized_decimal=True,
        source_enclosure_exact_mpf_tuples=[list(top),list(top)],
        positive_error_bound_atom_only=True,source_value_selected=False)
    return dict(name='O3_power_remaining_lower_finite_N_F_endpoint_error',
        coefficient_C0_C1=[_unit_row(),_unit_row()],log_scale_terms=[logatom],
        actual_function_binding=dict(chart='O3_power',native_coordinate='0',
            exact_signed_quantity='-R0^2*delta_F_lower',ordinary_Z_orders=[0,1],
            error_enclosure_only=True,source_budget_sample_sha256=source_digest,
            common_C0_C1_log_upper=True,lower_F_source_numerically_evaluated=False,
            no_aggregate_positive_bound_subtraction=True,Qy_or_viscous_term_added=False))


def materialize_Rc178_view(original,entries,independent_target,*,source_proof):
    require(source_proof.get('passed') is True and source_proof.get('same_source_Rc_F_pair_exact_zero') is True
        and source_proof.get('compact_upper_delta_F_C0_C1_exact_zero') is True,'Actual source F endpoint theorems required')
    require(len(original)==183 and len(entries)==178 and
        [v['parent_reduced_index'] for v in entries]==ACTIVE_INDICES,'Original183/Rc178 scope required')
    require(source_proof['independent_target_digest']==quiet_digest(independent_target),'Target was changed or reintroduced')
    restored=copy.deepcopy(entries)
    replaced=[i for i,v in enumerate(restored) if v['operation']==OPERATION]
    require(len(replaced)==1,'Exactly one Rc endpoint replacement required')
    row=restored[replaced[0]]
    require(row['parent_reduced_index']==176 and row['parent_term_digest']==quiet_digest(original[176]) and
        row['contracted_parent_179_digest']==quiet_digest(original[179]) and
        row['source_proof_digest']==digest(source_proof),'Original endpoint provenance or theorem changed')
    term=copy.deepcopy(row['term'])
    restored[replaced[0]]=dict(parent_reduced_index=176,parent_term_digest=quiet_digest(original[176]),
        name=original[176]['name'],operation='unchanged_parent_reference')
    restored.insert(ACTIVE_INDICES.index(180),dict(parent_reduced_index=179,
        parent_term_digest=quiet_digest(original[179]),name=original[179]['name'],operation='unchanged_parent_reference'))
    expanded=materialize_quiet179_view(original,restored,independent_target)
    at=INDICES.index(176);expanded[at]=term;del expanded[INDICES.index(179)]
    require(len(expanded)==178 and expanded[:6]==original[:6],'Six signed original core/lower terms changed')
    return expanded


def build_views(root,parent,profile,join,quiet,budget):
    require(all(v.get('passed') is True and v.get('source_executed') is False for v in (profile,join,quiet,budget)),
        'Passed source-only proof inputs required')
    cert=budget['compact_upper2_relative_C5_source_certificate']
    require(join['same_source_Rc_boundary_pair_exact_zero'] is True and
        join['deterministic_recovery_congruence']['F_ordinary_C1_equal'] is True and
        profile['compact_upper_delta_F_C0_C1_exact_zero'] is True and
        cert['passed'] is True and cert['source_only'] is True and cert['relative_mean_C5_exact_zero'] is True and
        cert['delta_F_zero_requires_actual_velocity_and_pressure_recovery_proof'] is True,
        'Actual Rc/upper source identity required before budget contraction')
    require(quiet['negative_remaining_lower_R0_squared_deltaF_C0_C1'] and
        budget['contraction_applied'] is False,'Remaining signed lower F source required')
    archive=parent['inherited_m15_packet'];path=Path(root)/'experiments/root_st073'/archive['report']
    require(sha(path)==archive['report_sha256'],'Original macro archive changed')
    raw=load(path);out=[]
    for old,saved,sample in zip(parent['quiet179_discrepancy_views'],raw['samples'],budget['samples']):
        require(old['Z']==saved['Z']==sample['Z'],'Different axial samples')
        original=saved['correlated_reduced_m15_C0_C1_terms'];target=saved['independent_target']
        require(quiet_digest(original)==old['parent_original_reduced183_digest'],'Archived reduced183 changed')
        prior=materialize_quiet179_view(original,old['entries'],target)
        proof=dict(passed=True,same_source_Rc_F_pair_exact_zero=True,
            compact_upper_delta_F_C0_C1_exact_zero=True,independent_target_digest=quiet_digest(target),
            Rc_pressure_join_sha256=digest(join),profile_join_sha256=digest(profile),
            remaining_lower_source_sha256=digest(quiet),budget_sample_sha256=digest(sample),
            exact_identity='[R^2*dF]_O3lower^Rc+[R^2*dF]_Rc^2Rc=-R0^2*dF(O3lower)',
            viscous_terms_already_contracted=True,viscous_external_term_untouched=True,
            all_volume_caps_retained=True,no_aggregate_positive_bound_subtraction=True)
        term=lower_endpoint_term(sample,digest(sample));entries=[]
        for entry in old['entries']:
            index=entry['parent_reduced_index']
            if index==179:continue
            if index==176:
                require(entry['operation']=='unchanged_parent_reference','O3 endpoint already replaced')
                entry=dict(entry,operation=OPERATION,contracted_parent_179_digest=quiet_digest(original[179]),
                    source_proof_digest=digest(proof),term=term)
            entries.append(copy.deepcopy(entry))
        expanded=materialize_Rc178_view(original,entries,target,source_proof=proof)
        require(all(expanded[j]==prior[INDICES.index(i)] for j,i in enumerate(ACTIVE_INDICES) if i!=176),
            'Unrelated active ledger term changed')
        out.append(dict(Z=old['Z'],parent_quiet179_view_sha256=digest(old),
            parent_original_reduced183_digest=quiet_digest(original),source_proof=proof,entries=entries,view_terms=178,
            expanded_Rc178_digest=digest(expanded),replaced_parent_indices=[176,179],
            archived_aggregate_endpoint_caps=[dict(index=i,term=original[i]) for i in (176,179)],
            all_other_177_terms_unchanged=True,independent_target_consumed_once=True,quiet_pair_consumed_once=True,
            signed_discrepancy_recipe='sum(Rc178 leaves); target and quiet pair already consumed',
            numeric_terminal_closure=False))
    require(len(out)==3,'Expected exactly three inherited diagnostic samples')
    return out
