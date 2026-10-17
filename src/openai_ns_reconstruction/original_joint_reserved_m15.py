"""Consume the common full/pulse C1 source at the selected full-support row.

Only reserved minus its independent negative-pulse target is contracted.
Every forward source category, including all finite-N errors, is retained.
The resulting interval is an enclosure, never a selected function value.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_hydro_source_point import ROOT,DIRECTORY,sha,digest,load,exact_iv
from .original_joint_m15_source import _ast_guard

PARENT='lei_ren_part1_paper_compliant_current_original_first_switch_terminal_m15_numeric.json'
PARENT_SHA256='1a3c4759f475b9eb615352ead770baa2b874f2e03485d2f4828eb0744ff36535'
FULL='lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge.json.gz'
FULL_SOURCE='lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge.py'
ORIGINAL_INDICES=[i for i in range(183) if i not in (38,39,40)]
JOINT_OPERATION='same_source_reserved_minus_target_contraction'


def record(ctx,value):
    return {'lower':ctx.nstr(ctx.make_mpf(value._mpi_[0]),45),
            'upper':ctx.nstr(ctx.make_mpf(value._mpi_[1]),45),
            'source_enclosure_exact_mpf_tuples':[list(v) for v in value._mpi_]}


def negate_component_sum(ctx,iv,components,*,source_bound):
    """Enclose the already contracted source, preserving every component box."""
    if source_bound is not True:raise ValueError('Exact common C1 source and full-support response required')
    if list(components)!=['inner','lower','preRp'] or any(len(v)!=2 for v in components.values()):
        raise ValueError('Exactly the three nonpulse ordinary C1 components required')
    return [record(ctx,-sum((exact_iv(ctx,iv,components[name][k])
        for name in ('inner','lower','preRp')),iv.mpf(0))) for k in (0,1)]


def contract_reserved_target(entries,joint,*,source_bound):
    if source_bound is not True:raise ValueError('Bound exact full-support response required')
    if len(entries)!=180 or [v['parent_reduced_index'] for v in entries]!=ORIGINAL_INDICES:
        raise ValueError('Accepted180 source view required')
    if sum(v['name']=='reserved_n1_own_swirl_once' for v in entries)!=1:
        raise ValueError('Exactly one original reserved term required')
    old=entries[-1]
    if old['parent_reduced_index']!=182 or old['name']!='reserved_n1_own_swirl_once':
        raise ValueError('Original reserved source slot182 required')
    if joint['parent_reduced_index']!=182 or joint['parent_term_digest']!=old['parent_term_digest']:
        raise ValueError('Joint source provenance changed')
    return copy.deepcopy(entries[:-1])+[copy.deepcopy(joint)]


def materialize_joint_discrepancy_view(original_terms,entries,independent_target):
    """Expand a discrepancy view: its external target is already consumed."""
    if len(original_terms)!=183 or len(entries)!=180 or [v['parent_reduced_index'] for v in entries]!=ORIGINAL_INDICES:
        raise ValueError('Original183 and joint180 source scope required')
    out=[];joint_count=0
    for row in entries:
        old=original_terms[row['parent_reduced_index']]
        if digest(old)!=row['parent_term_digest']:raise ValueError('Archived source term changed')
        if row['operation']=='unchanged_parent_reference':out.append(copy.deepcopy(old))
        elif row['operation']=='same_source_signed_leaf_replacement':out.append(copy.deepcopy(row['term']))
        elif row['operation']==JOINT_OPERATION:
            if row['parent_reduced_index']!=182 or row['independent_target_digest']!=digest(independent_target):
                raise ValueError('Different reserved/target source pair')
            if independent_target['actual_function_binding']['sign']!=-1:raise ValueError('Negative complete pulse target required')
            out.append(copy.deepcopy(row['term']));joint_count+=1
        else:raise ValueError('Unknown source-view operation')
    if joint_count!=1:raise ValueError('The independent target must be consumed exactly once')
    return out


class OriginalJointReservedM15:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Accepted terminal180 parent required')
        self.parent=load(self.directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed saved parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        self.verify_pins()
        self.ledger=load(self.directory/self.parent['inherited_m15_packet']['report'])
        self.full=load(self.directory/FULL)
        if self.full['source_family']!=self.parent['source_family']:raise ValueError('Different full pressure source family')
        self.ctx=mp.mp.clone();self.ctx.dps=240;self.iv=MPIntervalContext();self.iv.dps=280

    def verify_pins(self):
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed source pin: '+name)

    def cached_component_source(self):
        path=self.directory/FULL_SOURCE
        source={
            'sector_fixed_unit':_ast_guard(path,ROOT,'sector_C1_in_fixed_pressure_unit',
                statements=('bound=log-unit_log+c.ln(c.mpf(magnitude))',
                    'values[order]+=value'),
                text=('original_source_bases_are_fixed_within_Taylor_jets=True',
                    'source_ordinary_Z1_already_contains_F0_amplitude_derivatives=True',
                    'positive_unit_is_Z_independent=True','nonzero_sector_not_selected_zero=True')),
            'same_graph_components':_ast_guard(path,ROOT,'__init__',cls='CurrentOriginalN1FullFiveNumericBridge',
                statements=('self.functions=self.full.moments',),
                text=('o.g is not self.g','same positive reference*A*Rv power, Z independent')),
            'component_assembly':_ast_guard(path,ROOT,'evaluate',cls='CurrentOriginalN1FullFiveNumericBridge',
                statements=(
                    'inner_part=sector_C1_in_fixed_pressure_unit(c,sectors[key],self.unit_logs[key])',
                    "lower_part=sector_C1_in_fixed_pressure_unit(c,views['lower']['true_lower_pressure_C1_enclosures'][key],self.unit_logs[key])",
                    "other=[interval(c,v) for v in views['preRp']['normalized_true_preRp_C1_enclosures'][key]]",
                    "parts[key]=dict(inner=inner_part,lower=lower_part,preRp=other,pulse=last)"),
                text=('nonzero_source_errors_preserved=True','signed_m15_not_negated_again=True'))}
        if not all(v['passed'] is True for v in source.values()):raise ValueError('Cached normalized C1 component source binding failed')
        s=self.full['samples']['whole_Z'];parts=s['pressure_component_enclosures']['m15']
        components={name:copy.deepcopy(parts[name]['C1_enclosures'] if name!='preRp' else parts[name])
                    for name in ('inner','lower','preRp')}
        for name in ('inner','lower'):
            if not all(parts[name][k] is True for k in ('original_source_bases_are_fixed_within_Taylor_jets',
                'source_ordinary_Z1_already_contains_F0_amplitude_derivatives','positive_unit_is_Z_independent','nonzero_sector_not_selected_zero')):
                raise ValueError('Cached physical component errors/derivatives omitted')
        if self.full['source_identity']['pressure_unit_Z1_exact_zero'] is not True:
            raise ValueError('Fixed physical unit source theorem required')
        return {'passed':True,'AST':source,'report':FULL,'report_sha256':self.input_hashes[FULL],
            'source':FULL_SOURCE,'source_sha256':self.input_hashes[FULL_SOURCE],
            'Z_domain':copy.deepcopy(s['Z']),'ordinary_Z_orders':[0,1],
            'normalized_components':components,'normalized_components_digest':digest(components),
            'component_bounds_are_function_values':False,'component_bounds_or_midpoints_selected':False,
            'source_errors_and_amplitude_derivatives_preserved':True,'unit_Z_derivative_exact_zero':True}

    def build(self):
        from .original_reserved_full_support_source import build_reserved_full_support_source
        support=build_reserved_full_support_source(self.parent,root=ROOT,directory=self.directory)
        if support['passed'] is not True:raise ValueError('Full-support source binding failed: '+str(support.get('gaps')))
        for base,key,group in ((self.directory,'input_hashes',self.input_hashes),(ROOT,'implementation_hashes',self.implementation_hashes)):
            for name,expected in support[key].items():
                if group.get(name,expected)!=expected or sha(base/name)!=expected:raise ValueError('Conflicting support pin: '+name)
                group[name]=expected
        cached=self.cached_component_source()
        coeff=negate_component_sum(self.ctx,self.iv,cached['normalized_components'],source_bound=True)
        source=self.parent['joint_m15_pressure_source']
        joint_source={'scope':'selected full-support reserved response minus the independent negative complete-pulse m15 target only',
            'identity':'reserved-(-pulse)=-full+pulse=-(inner+lower+preRp)',
            'source_C1_pairs':copy.deepcopy(source['source_function_bindings']),
            'common_normalization':copy.deepcopy(source['normalized_expression']['shared_unit']),
            'ordinary_Z_orders':[0,1],'selected_compact_x':'2','exact_response':[0,0,1],
            'support_proof_digest':digest(support),'cached_component_source_digest':digest(cached),
            'source_identity_before_boxing':True,'independent_target_consumed_once':True,
            'coefficient_C0_C1':coeff,'enclosure_domain':copy.deepcopy(cached['Z_domain']),
            'enclosure_not_function_value':True,'source_centers_resolved':False,
            'physical_scale_materialized':False,'accepted_ancestor_or_quadrature_replayed':False}
        views=[]
        for before in self.parent['signed_m15_query_views']:
            target=before['independent_target']
            sample=next(v for v in self.ledger['samples'] if v['Z']==before['Z'])
            original=sample['correlated_reduced_m15_C0_C1_terms'];old=before['entries'][-1]
            term={'name':'reserved_minus_negative_pulse_joint_m15','coefficient_C0_C1':copy.deepcopy(coeff),
                'log_scale_terms':copy.deepcopy(target['log_scale_terms']),
                'actual_function_binding':copy.deepcopy(joint_source),
                'approximation_is_exact_source_value':False,'source_enclosure_is_not_a_selected_value':True}
            joint={'parent_reduced_index':182,'parent_term_digest':old['parent_term_digest'],
                'name':term['name'],'operation':JOINT_OPERATION,'term':term,
                'independent_target_digest':before['independent_target_digest'],
                'reserved_and_target_provenance':{'reserved':copy.deepcopy(old),
                    'independent_target':copy.deepcopy(target),'source_identity_proof_digest':digest(support)},
                'source_contraction_digest':digest(joint_source)}
            entries=contract_reserved_target(before['entries'],joint,source_bound=True)
            expanded=materialize_joint_discrepancy_view(original,entries,target)
            view=copy.deepcopy(before)
            view.update(entries=entries,view_semantics='physical discrepancy; independent target consumed inside the one joint source leaf',
                external_target_subtraction_required=False,independent_target_consumed_once=True,
                parent_forward180_entries_digest=digest(before['entries']),
                unchanged_parent_entries=179,unchanged_original_references=175,
                joint_reserved_minus_target_leaves=1,forward_six_inner_lower_cancelled=False,
                expanded_joint180_digest=digest(expanded),joint_source_digest=digest(joint_source),
                signed_discrepancy_recipe='sum(joint180 physical leaves); do not subtract independent target again',
                original_forward180_expansion_digest=before['expanded_signed180_digest'])
            # This old expansion describes the parent forward sum. It is kept
            # under an explicit provenance key, not mislabeled as this view.
            view.pop('expanded_signed180_digest')
            views.append(view)
        return {'source_family':copy.deepcopy(self.parent['source_family']),'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'reserved_full_support_source':support,'cached_component_source':cached,'joint_reserved_minus_target_source':joint_source,
            'joint_discrepancy_views':views,'inherited_m15_packet':copy.deepcopy(self.parent['inherited_m15_packet']),
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_forward180_views_digest':digest(self.parent['signed_m15_query_views']),
            'inherited_leading_packets_digest':digest(self.parent['macro_first_switch_leading_viscous_packets']),
            'contract':{'typed_full_support_joint_input_consumed':True,'joint_source_enclosure_computed':True,
                'shared_pulse_cancelled_as_function_before_boxing':True,'reserved_and_target_consumed_once':True,
                'other179_parent_entries_unchanged':True,'raw213_reduced183_archives_unchanged':True,
                'forward_six_inner_lower_and_all_finite_N_caps_retained':True,'exact_N_and_physical_units_retained':True,
                'accepted_ancestor_or_fixture_replayed':False,'physical_width_or_amplitude_materialized':False,
                'caps_midpoints_or_overlap_selected_as_values':False,'independent_target_subtracted_twice':False},
            'limits':{'signed_source_function_centers_resolved':False,'partial_response_input_correlations_resolved':False,
                'forward_canonical_component_function_equality':False,'numeric_Rp_or_global_five_closure':False,
                'exact_heat_energy_stress_flat':False,'true_temporal_n_recursion':False,
                'mean_oscillatory_correction':False,'independent_corrected_Cartesian_NS':False,'full_goal_complete':False}}
