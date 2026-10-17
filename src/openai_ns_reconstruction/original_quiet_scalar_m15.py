"""Pair the forward quiet increment with its canonical source exactly once.

The scalar identity is proved from the two source integrals. Independent
saved quadrature enclosures remain numerical diagnostics, never its proof.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_hydro_source_point import ROOT,DIRECTORY,sha,digest,load,exact_iv
from .original_joint_m15_source import _ast_guard
from .original_quiet_scalar_kernel import quiet_scalar_difference

PARENT='lei_ren_part1_paper_compliant_current_original_joint_reserved_m15_numeric.json'
PARENT_SHA256='a46413fd755aa96c385aa7b81e1919fe0d22838868b8b5a7d4348708dd16dad6'
AXIAL='lei_ren_part1_paper_compliant_current_original_axial_m15_numeric.json'
QUIET='lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source.json.gz'
PRERP='lei_ren_part1_paper_compliant_current_original_n1_preRp_numeric_bridge.json.gz'
PRERP_SOURCE='lei_ren_part1_paper_compliant_current_original_n1_preRp_numeric_bridge.py'
FULL='lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge.json.gz'
INDICES=[i for i in range(183) if i not in (38,39,40,181)]
REMAINDER_OPERATION='same_source_quiet_pair_contraction_remainder'


def record(ctx,value):
    return {'lower':ctx.nstr(ctx.make_mpf(value._mpi_[0]),45),
            'upper':ctx.nstr(ctx.make_mpf(value._mpi_[1]),45),
            'source_enclosure_exact_mpf_tuples':[list(v) for v in value._mpi_]}


def contract_quiet_pair(entries,remainder,*,source_equality_proved):
    if source_equality_proved is not True:raise ValueError('Same integral and physical K source theorem required')
    expected=[i for i in range(183) if i not in (38,39,40)]
    if len(entries)!=180 or [v['parent_reduced_index'] for v in entries]!=expected:
        raise ValueError('Accepted joint180 discrepancy scope required')
    if entries[-2]['parent_reduced_index']!=181 or entries[-2]['name']!='quiet_exact_original_mean_square_integral':
        raise ValueError('Original quiet increment required')
    if entries[-1]['operation']!='same_source_reserved_minus_target_contraction':
        raise ValueError('Target must already be consumed by the accepted joint source')
    if remainder['parent_term_digest']!=entries[-1]['parent_term_digest']:
        raise ValueError('Original reserved provenance changed')
    return copy.deepcopy(entries[:-2])+[copy.deepcopy(remainder)]


def materialize_quiet179_view(original,entries,independent_target):
    if len(original)!=183 or len(entries)!=179 or [v['parent_reduced_index'] for v in entries]!=INDICES:
        raise ValueError('Original183 and quiet179 scope required')
    out=[];joined=0
    for row in entries:
        old=original[row['parent_reduced_index']]
        if digest(old)!=row['parent_term_digest']:raise ValueError('Archived source changed')
        if row['operation']=='unchanged_parent_reference':out.append(copy.deepcopy(old))
        elif row['operation']=='same_source_signed_leaf_replacement':out.append(copy.deepcopy(row['term']))
        elif row['operation']==REMAINDER_OPERATION:
            if row['parent_reduced_index']!=182 or row['independent_target_digest']!=digest(independent_target):
                raise ValueError('Changed source/target group')
            out.append(copy.deepcopy(row['term']));joined+=1
        else:raise ValueError('Unknown source operation')
    if joined!=1:raise ValueError('Exactly one source/target remainder required')
    return out


class OriginalQuietScalarM15:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Accepted joint source parent required')
        self.parent=load(self.directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed joint source parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed saved source pin: '+name)
        self.ledger=load(self.directory/self.parent['inherited_m15_packet']['report'])
        self.ctx=mp.mp.clone();self.ctx.dps=240;self.iv=MPIntervalContext();self.iv.dps=280

    def independent_diagnostic(self,proof):
        c,iv=self.ctx,self.iv
        axial=load(self.directory/AXIAL)
        ends=[v for v in axial['pressure_row_samples'] if v['native_coordinate']=='1']
        if not ends:raise ValueError('Saved axial terminal mean source missing')
        mean=ends[0]['leading_scalars']['mean']
        if any(v['leading_scalars']['mean']!=mean for v in ends):
            raise ValueError('One Z-independent saved axial terminal mean required')
        raw=load(self.directory/QUIET);q=raw['actual_quiet_source']
        logF=q['kernel']['log_F'];length=q['original_quiet_length']
        a=4*exact_iv(c,iv,mean);b=4*iv.exp(exact_iv(c,iv,logF));L=exact_iv(c,iv,length)
        if a.a<=0 or b.a<=0:raise ValueError('Positive original scalar means required')
        uncorrelated,receipt=quiet_scalar_difference(iv,a,b,L,source_bound=True)
        correlated,identity_receipt=quiet_scalar_difference(iv,a,b,L,source_bound=True,
            gap=iv.mpf(0),gap_source_bound=proof['source_mean_equality']['passed'])
        if any(v._mpi_!=iv.mpf(0)._mpi_ for v in correlated):raise ValueError('Exact scalar identity did not contract')
        return {'ordinary_Z_orders':[0,1],'source_mean_summary':copy.deepcopy(mean),
            'canonical_log_F':copy.deepcopy(logF),'physical_log_length':copy.deepcopy(length),
            'forward_K_over_Raxial_C0':record(c,a),'canonical_K_over_Raxial_C0':record(c,b),
            'independent_saved_enclosure_difference_C0_C1':[record(c,v) for v in uncorrelated],
            'same_source_correlated_difference_C0_C1':[record(c,v) for v in correlated],
            'independent_kernel_receipt':receipt,'source_identity_kernel_receipt':identity_receipt,
            'common_physical_log_scale':'2log(R_axial_end), copied from original quiet term',
            'normalized_scalar_exponential_only':True,'physical_radius_width_or_amplitude_materialized':False,
            'mean_approximation_or_bound_selected_as_actual_value':False,
            'interval_overlap_used_as_source_identity':False,'same_source_integral_theorem_used':True,
            'accepted_quadrature_or_owner_replayed':False}

    def remaining_components(self):
        full=load(self.directory/FULL)['samples']['whole_Z']
        pre=load(self.directory/PRERP)['samples']['whole_Z']
        rows=pre['cell_piece_enclosures'];c,iv=self.ctx,self.iv
        if not rows or any(set(v['piece_enclosures']['m15'])!={'active','compact','quiet'} for v in rows):
            raise ValueError('All canonical preRp pieces required')
        comp={name:copy.deepcopy(full['pressure_component_enclosures']['m15'][name]['C1_enclosures']) for name in ('inner','lower')}
        for name in ('active','compact'):
            values=[[exact_iv(c,iv,v['piece_enclosures']['m15'][name][k]) for v in rows] for k in (0,1)]
            comp[name]=[record(c,iv.mpf([min(v.a for v in values[k]),max(v.b for v in values[k])])) for k in (0,1)]
        source=_ast_guard(self.directory/PRERP_SOURCE,ROOT,'evaluate',cls='CurrentOriginalN1PreRpNumericBridge',
            statements=("pieces[key]=dict(active=[negative,negative],compact=[negative,negative],quiet=[quiet,c.mpf(0)])",),
            text=('ordinary_Z0_Z1_enclosures_not_defining_function_values=True','inner_lower_and_pulse_pressure_not_added=True'))
        if source['passed'] is not True:raise ValueError('Canonical piece bound source changed')
        unit_keys=('source_logRc_formula','source_logRp_formula','source_Rv_formula','source_A',
                   'source_logE_formula','normalization','m15')
        if any(full['pressure_normalization'][k]!=pre['normalization'][k] for k in unit_keys):
            raise ValueError('Nonquiet canonical units changed')
        coeff=[record(c,-sum((exact_iv(c,iv,v[k]) for v in comp.values()),iv.mpf(0))) for k in (0,1)]
        return {'passed':True,'AST':{'preRp_piece_enclosures':source},'ordinary_Z_orders':[0,1],
            'components':comp,'component_names':['inner','lower','active','compact'],
            'coefficient_C0_C1':coeff,'enclosure_domain':copy.deepcopy(full['Z']),
            'canonical_physical_C1_pairs':{'inner':[53161,53162],'lower':[36592,36593],
                'active':[30707,30708],'compact':[31155,31156]},
            'source_rows_are_enclosures_not_function_values':True,'source_errors_retained':True,
            'cached_full_report':FULL,'cached_preRp_report':PRERP}

    def build(self):
        from .original_quiet_scalar_source import build_quiet_scalar_source
        proof=build_quiet_scalar_source(self.parent,root=ROOT,directory=self.directory)
        if proof['passed'] is not True or proof['source_mean_equality']['passed'] is not True:
            raise ValueError('Quiet same-source scalar identity required: '+str(proof.get('gaps')))
        for base,key,group in ((self.directory,'input_hashes',self.input_hashes),(ROOT,'implementation_hashes',self.implementation_hashes)):
            for name,expected in proof[key].items():
                if group.get(name,expected)!=expected or sha(base/name)!=expected:raise ValueError('Conflicting quiet source pin: '+name)
                group[name]=expected
        diagnostic=self.independent_diagnostic(proof);remaining=self.remaining_components();views=[]
        binding={'identity':'forward_quiet-canonical_quiet=0; remaining joint=-(inner+lower+active+compact)',
            'canonical_quiet_pair':[31215,0],'ordinary_Z_orders':[0,1],
            'source_mean_equality_proof_digest':digest(proof),
            'canonical_nonquiet_source_digest':digest(remaining),
            'physical_unit':'unchanged independent-target reference*A*Rv^2, U_Z=0',
            'source_correlated_before_boxing':True,'signed_source_centers_resolved':False,
            'independent_target_already_consumed_once':True,'numerical_enclosure_not_function_value':True}
        for before in self.parent['joint_discrepancy_views']:
            sample=next(v for v in self.ledger['samples'] if v['Z']==before['Z'])
            original=sample['correlated_reduced_m15_C0_C1_terms'];quiet=original[181]
            old=before['entries'][-1];target=before['independent_target']
            term={'name':'reserved_target_after_exact_quiet_pair','coefficient_C0_C1':copy.deepcopy(remaining['coefficient_C0_C1']),
                'log_scale_terms':copy.deepcopy(target['log_scale_terms']),'actual_function_binding':copy.deepcopy(binding),
                'approximation_is_exact_source_value':False,'source_enclosure_is_not_a_selected_value':True}
            remainder={'parent_reduced_index':182,'parent_term_digest':old['parent_term_digest'],
                'name':term['name'],'operation':REMAINDER_OPERATION,'term':term,
                'independent_target_digest':before['independent_target_digest'],
                'parent_joint_source_entry_digest':digest(old),'source_mean_equality_proof_digest':digest(proof)}
            entries=contract_quiet_pair(before['entries'],remainder,source_equality_proved=True)
            expanded=materialize_quiet179_view(original,entries,target)
            view=copy.deepcopy(before)
            view.update(entries=entries,view_terms=179,unchanged_parent_entries=178,unchanged_original_references=174,
                view_semantics='discrepancy; target already consumed; actual quiet source pair cancelled',
                parent_joint180_entries_digest=digest(before['entries']),parent_joint180_expansion_digest=before['expanded_joint180_digest'],
                expanded_quiet179_digest=digest(expanded),quiet_scalar_diagnostic_digest=digest(diagnostic),
                exact_zero_quiet_source_provenance={'parent_entry':copy.deepcopy(before['entries'][-2]),
                    'original_quiet_term':copy.deepcopy(quiet),'canonical_quiet_pair':[31215,0],
                    'source_mean_equality_proof_digest':digest(proof),'actual_source_zero':True,
                    'old_approximation_enclosure_retained_as_provenance':True,'old_bound_selected_as_value':False},
                quiet_pair_consumed_once=True,signed_discrepancy_recipe='sum(quiet179 physical leaves); target and exact quiet pair already consumed')
            view['inherited_reserved_target_joint_source_digest']=view.pop('joint_source_digest')
            view['nonquiet_remainder_source_digest']=digest(remaining)
            view.pop('expanded_joint180_digest');views.append(view)
        return {'source_family':copy.deepcopy(self.parent['source_family']),'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'quiet_scalar_source':proof,'independent_scalar_diagnostic':diagnostic,'canonical_nonquiet_remainder':remaining,
            'quiet179_discrepancy_views':views,'inherited_m15_packet':copy.deepcopy(self.parent['inherited_m15_packet']),
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_joint180_views_digest':digest(self.parent['joint_discrepancy_views']),
            'inherited_leading_packets_digest':self.parent['inherited_leading_packets_digest'],
            'contract':{'same_source_integral_and_physical_K_identity':True,'quiet_pair_consumed_once':True,
                'target_already_consumed_once':True,'other178_parent_entries_unchanged':True,
                'original_raw213_reduced183_six_signed_inner_lower_and_N_unchanged':True,
                'all_finite_N_caps_and_pressure_tail_retained':True,
                'saved_independent_quadrature_diagnostic_retained':True,
                'numerical_overlap_midpoint_or_cap_used_as_source_value':False,
                'accepted_ancestor_or_quadrature_replayed':False,'physical_radius_width_amplitude_materialized':False},
            'limits':{'remaining_signed_source_centers_resolved':False,'forward_core_collar_active_compact_equality':False,
                'numeric_Rp_or_global_five_closure':False,'exact_heat_energy_stress_flat':False,'true_temporal_n_recursion':False,
                'mean_oscillatory_correction':False,'independent_corrected_Cartesian_NS':False,'full_goal_complete':False}}
