"""Source-bound first-switch terminal E contraction over saved m15 data.

Only the unmodified Uy endpoint vanishes. The own mean, H=R*Qy, finite-N
corrections and all other moment categories remain present.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_abs,mpf_cmp

from .original_hydro_source_point import ROOT,DIRECTORY,sha,digest,load,exact_iv

PARENT='lei_ren_part1_paper_compliant_current_original_macro_coupled_m15_numeric.json'
PARENT_SHA256='b63476443d5de5c1fca24ed5979ccb8a71273c42ff7a0b9b9094d8f80782f177'
ORIGINAL_INDICES=[i for i in range(183) if i not in (38,39,40)]
LEADING_SLOTS=(19,20,41,42)


def contract_endpoint_entries(entries,*,source_zero_proved):
    """Contract only the two original endpoint categories after proof."""
    if source_zero_proved is not True:raise ValueError('Exact truth-source endpoint theorem required')
    if len(entries)!=182 or [v['parent_reduced_index'] for v in entries]!=[i for i in range(183) if i!=39]:
        raise ValueError('Accepted signed182 entries required')
    kept=[];removed=[]
    for row in entries:
        index=row['parent_reduced_index']
        if index in (38,40):
            raw=43 if index==38 else 45
            removed.append(dict(copy.deepcopy(row),raw_index=raw,active=False,provenance_only=True,
                exact_actual_source_function_zero=True,old_enclosure_selected_as_value=False))
        else:kept.append(copy.deepcopy(row))
    return kept,removed


def materialize_terminal_view(original_terms,entries):
    """Expand original reduced183 references with the terminal contractions."""
    if len(original_terms)!=183 or len(entries)!=180:raise ValueError('Expected original183 and terminal180')
    if [v['parent_reduced_index'] for v in entries]!=ORIGINAL_INDICES:raise ValueError('Wrong original endpoint scope')
    out=[]
    for row in entries:
        old=original_terms[row['parent_reduced_index']]
        if digest(old)!=row['parent_term_digest']:raise ValueError('Changed archived source term')
        if row['operation']=='unchanged_parent_reference':out.append(copy.deepcopy(old))
        elif row['operation']=='same_source_signed_leaf_replacement':out.append(copy.deepcopy(row['term']))
        else:raise ValueError('Unknown source view operation')
    return out


class OriginalFirstSwitchTerminalM15:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Accepted coupled source packet required')
        self.parent=load(self.directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed saved parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted source pin: '+name)
        self.ledger=load(self.directory/self.parent['inherited_m15_packet']['report'])

    def leading_query_bounds(self,components,target):
        """Bound physical ordinary rows in target units without exponentiation.

        This is a bound on |D_Z^k leading|/unit(Z), not a derivative of a
        normalized ratio or a selected function value. Every physical scale
        remains logarithmic, including the original width and amplitude.
        """
        c=mp.mp.clone();c.dps=180;iv=MPIntervalContext();iv.dps=220
        unit=sum((exact_iv(c,iv,v) for v in target['log_scale_terms']),iv.mpf(0))
        rows=[]
        for order in (0,1):
            caps=[]
            for component in components:
                term=component['term'];pair=term['coefficient_C0_C1'][order]
                ends=[mpf_abs(tuple(v)) for v in pair['source_enclosure_exact_mpf_tuples']]
                top=ends[0] if mpf_cmp(ends[0],ends[1])>=0 else ends[1]
                if top[1]==0:continue
                logscale=sum((exact_iv(c,iv,v) for v in term['log_scale_terms']),iv.mpf(0))
                bound=iv.log(iv.mpf(c.make_mpf(top)))+logscale-unit
                caps.append(bound._mpi_[1])
            top=max(caps,key=lambda v:c.make_mpf(v)) if caps else None
            bound=iv.mpf(c.make_mpf(top))+iv.log(len(caps)) if caps else None
            rows.append({'ordinary_Z_order':order,'nonzero_component_enclosures':len(caps),
                'normalized_absolute_upper_log_exact_mpf_tuple':list(bound._mpi_[1]) if bound is not None else None,
                'normalized_absolute_upper_log_approximation':c.nstr(c.make_mpf(bound._mpi_[1]),35) if bound is not None else '-inf',
                'source_value_selected':False,'bound_only_not_function_value':True})
        return {'rows':rows,'normalization':'|D_Z^k leading physical m15| / independent target reference*A*Rv^2 at this query',
            'normalization_log_atoms_digest':digest(target['log_scale_terms']),
            'sum_bound_method':'log(number of nonzero enclosures)+max(component log upper)',
            'width_or_amplitude_exponentiated':False,'NS_or_terminal_validation_gate':False,
            'independent_target_coefficient_C0_C1':copy.deepcopy(target['coefficient_C0_C1']),
            'independent_target_coefficient_is_not_a_selected_value':True,
            'three_query_bounds_do_not_prove_global_axial_closure':True}

    def build(self):
        from .original_first_switch_zero_source import build_first_switch_zero_source
        proof=build_first_switch_zero_source(self.parent,root=ROOT,directory=self.directory)
        if proof['passed'] is not True:raise ValueError('First-switch truth-source zero theorem failed: '+str(proof.get('gaps')))
        from .original_joint_m15_source import build_joint_m15_source
        joint=build_joint_m15_source(self.parent,root=ROOT,directory=self.directory)
        if joint['passed'] is not True:raise ValueError('Common C1 pressure source bundle failed: '+str(joint.get('gaps')))
        for base,key,group in ((self.directory,'input_hashes',self.input_hashes),
                               (ROOT,'implementation_hashes',self.implementation_hashes)):
            for name,expected in joint[key].items():
                if group.get(name,expected)!=expected or sha(base/name)!=expected:
                    raise ValueError('Conflicting joint-source pin: '+name)
                group[name]=expected
        views=[];packets=[]
        for before in self.parent['signed_m15_query_views']:
            sample=next(v for v in self.ledger['samples'] if v['Z']==before['Z'])
            original=sample['correlated_reduced_m15_C0_C1_terms']
            entries,removed=contract_endpoint_entries(before['entries'],source_zero_proved=proof['passed'])
            for v in removed:v['actual_source_zero_proof_digest']=digest(proof)
            expanded=materialize_terminal_view(original,entries)
            by_index={v['parent_reduced_index']:v for v in entries}
            packet={'Z':before['Z'],'scope':'unmodified macro plus complete first-switch leading viscous category only',
                'component_parent_slots':list(LEADING_SLOTS),
                'components':[copy.deepcopy(by_index[i]) for i in LEADING_SLOTS],
                'identity':'Hnative-Hfirst1=Enative+I_macro+I_first; Efirst1=0',
                'first_switch_source_integral':'I_first=100*epsilon*hbS^2*C(sumBj)*M0+100*epsilon*hbS^3*C(R)',
                'physical_m15_units':'unchanged physical moment ledger; independent target reference*A*Rv^2',
                'ordinary_Z_orders':[0,1],'first1_H_or_absolute_mean_zero_claimed':False,
                'total_m15_or_terminal_closure_claimed':False,
                'independent_target_digest':before['independent_target_digest'],
                'independent_target_external_multiplier_in_discrepancy':-1,
                'full_discrepancy_partition':'leading_packet+sum(other176 original references)-independent_negative_full_pulse_m15_target'}
            packet['query_bound_comparison']=self.leading_query_bounds(packet['components'],before['independent_target'])
            packets.append(packet)
            views.append({'Z':before['Z'],'raw213_digest':before['raw213_digest'],
                'parent_original_reduced183_digest':before['parent_original_reduced183_digest'],
                'parent_signed182_view_digest':before['expanded_signed182_digest'],
                'entries':entries,'view_terms':180,'unchanged_original_references':176,
                'four_signed_replacement_slots':list(LEADING_SLOTS),
                'contracted_parent_slots':[38,39,40],'new_exact_zero_parent_slots':[38,40],
                'inherited_contracted_raw44_provenance':copy.deepcopy(before['contracted_source_endpoint_provenance']),
                'exact_zero_endpoint_provenance':removed,'expanded_signed180_digest':digest(expanded),
                'leading_viscous_packet_digest':digest(packet),
                'independent_target':copy.deepcopy(before['independent_target']),
                'independent_target_digest':before['independent_target_digest'],
                'selected_N_evidence_digest':before['selected_N_evidence_digest'],
                'six_inner_lower_digest':before['six_inner_lower_digest'],
                'largest_unchanged_caps':copy.deepcopy(before['largest_unchanged_caps']),
                'signed_discrepancy_recipe':'sum(signed180 physical leaves)-unchanged independent negative complete-pulse m15 target',
                'global_terminal_closure':False})
        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'first_switch_truth_source_zero':proof,'signed_m15_query_views':views,
            'joint_m15_pressure_source':joint,'joint_pressure_source_used_as_total_reserved_value':False,
            'macro_first_switch_leading_viscous_packets':packets,
            'inherited_m15_packet':self.parent['inherited_m15_packet'],
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_signed182_views_digest':digest(self.parent['signed_m15_query_views']),
            'inherited_macro_weighted_source_rows_digest':self.parent['macro_weighted_source_rows_digest'],
            'inherited_native_macro_seam_digest':digest(self.parent['native_macro_seam']),
            'inherited_signed_ledger_binding_digest':digest(self.parent['signed_ledger_binding']),
            'contract':{'accepted_ancestor_or_fixture_replayed':False,'accepted_source_integrals_or_defect_rows_recomputed':False,
                'raw213_reduced183_archives_unchanged':True,'prior_four_signed_leaves_unchanged':True,
                'six_inner_lower_reserved_once_target_and_N_unchanged':True,
                'positive_caps_midpoints_rounded_zero_or_interval_overlap_used_as_values':False,
                'hbB_hbS_or_absolute_swirl_amplitude_materialized':False,
                'first1_H_mean_or_finite_N_fields_zeroed':False,'all_other176_original_references_retained':True},
            'limits':{'numeric_global_axial_function_primitive':False,'numeric_Rp_terminal_join':False,
                'global_five_moment_common_N_closure':False,'exact_heat_energy_stress_flat':False,
                'true_temporal_n_recursion':False,'mean_oscillatory_correction':False,
                'independent_corrected_Cartesian_NS':False,'full_goal_complete':False}}
