"""Produce the source-zero terminal view, or inspect saved bytes only."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_first_switch_terminal_m15 import (
    PARENT,PARENT_SHA256,OriginalFirstSwitchTerminalM15,contract_endpoint_entries,
    materialize_terminal_view,sha,digest,load,LEADING_SLOTS,
)

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_first_switch_terminal_m15_numeric.json'
CODE=('src/openai_ns_reconstruction/original_first_switch_zero_source.py',
      'src/openai_ns_reconstruction/original_joint_m15_source.py',
      'src/openai_ns_reconstruction/original_first_switch_terminal_m15.py',
      'scripts/run_original_first_switch_terminal_m15.py')


def require(value,message):
    if not value:raise AssertionError(message)


def assert_unselected(value):
    if isinstance(value,dict):
        if 'approximation_is_exact_source_value' in value:require(value['approximation_is_exact_source_value'] is False,'Approximation selected')
        if 'source_enclosure_is_not_a_selected_value' in value:require(value['source_enclosure_is_not_a_selected_value'] is True,'Enclosure selected')
        for v in value.values():assert_unselected(v)
    elif isinstance(value,(list,tuple)):
        for v in value:assert_unselected(v)


def new_controls():
    """Exact finite algebra; no saved source or accepted fixture is run."""
    controls=[]
    def check(name,condition):
        require(condition,name);controls.append({'name':name,'passed':True})
    def endpoint(s,z):
        # A nonlinear axial source and a finite cutoff model test the zero
        # factor algebra. Actual original cutoff is bound by the AST theorem.
        sigma=3*s*s-2*s*s*s;k=1-sigma;delta=F(1,7)
        b=k*(2+z+z*z);bz=k*(1+2*z);den=1-delta*z*z
        return (2*z*b/den,2*(b+z*bz)/den+4*delta*z*z*b/(den*den))
    for z in (F(-1,3),F(0),F(1,3)):
        a=endpoint(F(1),z)
        check('exact_source_factor_Z0_'+str(z),a[0]==0)
        check('ordinary_Z1_and_actual_radius_factor_'+str(z),F(137,100)*a[1]==0)
    a=endpoint(F(1,2),F(0))
    check('point_E0_zero_does_not_prove_source_or_Z1_zero',a[0]==0 and a[1]!=0)
    check('near_endpoint_nonzero_is_not_flat_zero',endpoint(F(999,1000),F(1,3))[0]!=0)
    entries=[{'parent_reduced_index':i} for i in range(183) if i!=39]
    for bad in (False,1):
        try:contract_endpoint_entries(entries,source_zero_proved=bad)
        except ValueError:check('missing_exact_boolean_source_proof_'+str(bad),True)
        else:raise AssertionError('Unproved source zero accepted')
    try:contract_endpoint_entries(entries[:-1],source_zero_proved=True)
    except ValueError:check('wrong_original_source_scope_rejected',True)
    else:raise AssertionError('Wrong source scope accepted')
    # A shared nonconstant unit must carry its ordinary Z derivative.
    z=F(1,3);u=1+z*z;uz=2*z
    pre=(z+z*z+z*z*z,1+2*z+3*z*z)
    pulse=(2+3*z-z**4,3-4*z**3);full=tuple(a+b for a,b in zip(pre,pulse))
    normalize=lambda q:(q[0]/u,q[1]/u-q[0]*uz/(u*u))
    a,b,c=normalize(full),normalize(pulse),normalize(pre)
    joint=tuple(-x+y for x,y in zip(a,b))
    check('shared_C1_unit_cancellation_Z0',joint[0]==-c[0])
    check('shared_C1_unit_cancellation_Z1',joint[1]==-c[1])
    check('unit_Z_derivative_cannot_be_discarded',joint[1]!=-pre[1]/u)
    return {'passed':True,'control_count':len(controls),'controls':controls,
        'new_exact_algebra_only':True,'original_cutoff_model_executed':False,
        'accepted_owner_fixture_or_quadrature_replayed':False,
        'fixture_is_not_original_cutoff_source_proof':True,'shared_unit_C1_product_rule_included':True}


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Pinned coupled parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'],'Passed parent required');return p


def check_packet(r,p):
    require(r['source_family']==p['source_family'] and r['parent_report_sha256']==PARENT_SHA256,'Family/parent changed')
    proof=r['first_switch_truth_source_zero']
    require(proof['passed'] is True,'Actual truth-source endpoint theorem failed')
    require(all(v['passed'] is True for v in proof['AST'].values()),'Endpoint source AST failed')
    for key in ('sourceflagtruth','sigmaendpointbranch','ordinaryZsourcezeroorders','scopephase1','capnotselected','finiteN_meannotzeroed'):
        require(proof[key]['passed'] is True,'Actual source boundary failed: '+key)
    require(proof['ordinaryZsourcezeroorders']['orders']==[0,1] and proof['scopephase1']['phase']=='s=1','Wrong source phase/order')
    require(proof['sourceflagtruth']['truth_source_multiplied_before_ordinary_Z'] is True,'Only a numerical point is zero')
    require(proof['zero_scope']['actual_mean_zero_claim'] is False and proof['zero_scope']['H_field_zero_claim'] is False,'Nonzero mean/H zeroed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'Original witness changed')
    require(r['inherited_signed182_views_digest']==digest(p['signed_m15_query_views']),'Parent views changed')
    require(r['inherited_macro_weighted_source_rows_digest']==p['macro_weighted_source_rows_digest'],'Saved primitive changed')
    require(r['inherited_native_macro_seam_digest']==digest(p['native_macro_seam']) and r['inherited_signed_ledger_binding_digest']==digest(p['signed_ledger_binding']),'Accepted correlation changed')
    joint=r['joint_m15_pressure_source']
    require(joint['passed'] is True and all(v['passed'] is True for v in joint['AST'].values()),'Common C1 pressure source relation failed')
    require(joint['correlation_proved'] is True and joint['graph_identity']['passed'] is True,'Unproved producer partition')
    require(all(joint['graph_identity']['normalization_checks'].values()),'Common physical unit graph changed')
    require(joint['source_function_bindings']=={'inner':[53161,53162],'lower':[36592,36593],
        'preRp':[31221,31220],'pulse':[35556,35557],'full':[53175,53176]},'Wrong m15 function partition')
    require(joint['normalized_expression']['ordinary_Z_orders']==[0,1] and
        joint['normalized_expression']['shared_unit']['U_Z_exact_zero'] is True,'Ordinary C1 source scaling changed')
    require(joint['graph_identity']['free_bounded_atom_used'] is False,'Independent bounded atom substituted')
    require(all(joint[k] is False for k in ('numeric_value_available','numeric_joint_response_complete',
        'full_forward_partition_mapping_complete','forward_partition_consumer_bound','global_closure_claimed')),
        'Producer relation promoted to numerical consumer closure')
    require(all(v is False for v in joint['unresolved_consumer_center_angle_coupling'].values()),'Uncomputed consumer correlation accepted')
    require(all(joint['physical_scale_guards'][k] is True for k in ('Pstar_applied_once','physical_measure','signed_m15_once')),
        'Physical pressure measure/sign changed')
    require(joint['cached_bounds_refs']['bounds_are_function_values'] is False and
        joint['cached_bounds_refs']['source_enclosure_midpoint_used'] is False,'Box promoted to source value')
    for key in ('input_hashes','implementation_hashes'):
        require(all(r[key].get(name)==expected for name,expected in joint[key].items()),'Joint source pin missing from publication')
    forward_pin=joint['source_pins']['forward_native_pin_report']
    require(forward_pin['passed'] is True and forward_pin['expected_sha256']==
        '311c61f1bc38143498ba443acac266b7160001037ef4bc022e52d9faf52940f0','Native consumer lineage report is unanchored')
    require(r['joint_pressure_source_used_as_total_reserved_value'] is False,'Partial pressure input promoted to full response')
    ledger=load(DIRECTORY/p['inherited_m15_packet']['report'])
    require(len(r['signed_m15_query_views'])==len(r['macro_first_switch_leading_viscous_packets'])==3,'Three queries required')
    for view,before,packet in zip(r['signed_m15_query_views'],p['signed_m15_query_views'],r['macro_first_switch_leading_viscous_packets']):
        require(view['Z']==before['Z']==packet['Z'],'Query order changed')
        old=next(v for v in ledger['samples'] if v['Z']==view['Z'])
        terms=old['correlated_reduced_m15_C0_C1_terms'];raw=old['raw_m15_C0_C1_terms']
        require(len(terms)==183 and len(raw)==213 and digest(raw)==view['raw213_digest'],'Original archive changed')
        require(digest(terms)==view['parent_original_reduced183_digest'],'Original reduced183 changed')
        require(view['parent_signed182_view_digest']==before['expanded_signed182_digest'],'Parent signed view changed')
        require(view['view_terms']==len(view['entries'])==180 and view['unchanged_original_references']==176,'Terminal scope count changed')
        expanded=materialize_terminal_view(terms,view['entries'])
        require(digest(expanded)==view['expanded_signed180_digest'],'Terminal expansion changed')
        expected=[v for v in before['entries'] if v['parent_reduced_index'] not in (38,40)]
        require(view['entries']==expected,'A retained source entry changed')
        require(view['four_signed_replacement_slots']==list(LEADING_SLOTS),'Signed source slots changed')
        require(view['contracted_parent_slots']==[38,39,40] and view['new_exact_zero_parent_slots']==[38,40],'Wrong contraction scope')
        require(view['inherited_contracted_raw44_provenance']==before['contracted_source_endpoint_provenance'],'Previous raw44 provenance changed')
        removed=view['exact_zero_endpoint_provenance']
        require([(v['parent_reduced_index'],v['raw_index']) for v in removed]==[(38,43),(40,45)],'Wrong source endpoints retired')
        for v in removed:
            require(v['parent_term_digest']==digest(terms[v['parent_reduced_index']]),'Archived zero-source provenance changed')
            require(v['active'] is False and v['provenance_only'] and v['exact_actual_source_function_zero'],'Endpoint active or theorem missing')
            require(v['old_enclosure_selected_as_value'] is False and v['actual_source_zero_proof_digest']==digest(proof),'Cap promoted to source value')
        require(expanded[:6]==terms[:6] and view['six_inner_lower_digest']==digest(raw[:6]),'Six signed inner/lower changed')
        require(sum(v['name']=='reserved_n1_own_swirl_once' for v in expanded)==1,'Reserved count changed')
        require(view['independent_target']==old['independent_target'] and view['independent_target_digest']==digest(old['independent_target']),'Independent target changed')
        binding=view['independent_target']['actual_function_binding']
        require(binding['sign']==-1 and binding['physical_measure']=='-R*pressure_density*dR','P1/positive target substituted')
        require(view['selected_N_evidence_digest']==digest(old['selected_N_evidence']),'Common N changed')
        require(view['largest_unchanged_caps']==before['largest_unchanged_caps'] and view['global_terminal_closure'] is False,'Old limits changed')
        require(packet['component_parent_slots']==list(LEADING_SLOTS) and packet['ordinary_Z_orders']==[0,1],'Leading packet scope changed')
        byindex={v['parent_reduced_index']:v for v in view['entries']}
        require(packet['components']==[byindex[i] for i in LEADING_SLOTS],'Leading components changed')
        require(packet['first1_H_or_absolute_mean_zero_claimed'] is False and packet['total_m15_or_terminal_closure_claimed'] is False,'Mean or full moment zeroed')
        require(packet['independent_target_digest']==view['independent_target_digest'] and packet['independent_target_external_multiplier_in_discrepancy']==-1,'Discrepancy target changed')
        q=packet['query_bound_comparison']
        require(q['normalization_log_atoms_digest']==digest(view['independent_target']['log_scale_terms']),'Wrong physical query units')
        require(q['independent_target_coefficient_C0_C1']==view['independent_target']['coefficient_C0_C1'],'Independent target coefficient changed')
        require(q['width_or_amplitude_exponentiated'] is False and q['NS_or_terminal_validation_gate'] is False,'Bound promoted to validation')
        require([v['ordinary_Z_order'] for v in q['rows']]==[0,1] and all(v['bound_only_not_function_value'] and v['source_value_selected'] is False for v in q['rows']),'Bound used as function value')
        require(digest(packet)==view['leading_viscous_packet_digest'],'Leading packet digest changed')
        assert_unselected(view)
    true_keys=('raw213_reduced183_archives_unchanged','prior_four_signed_leaves_unchanged',
        'six_inner_lower_reserved_once_target_and_N_unchanged','all_other176_original_references_retained')
    false_keys=('accepted_ancestor_or_fixture_replayed','accepted_source_integrals_or_defect_rows_recomputed',
        'positive_caps_midpoints_rounded_zero_or_interval_overlap_used_as_values',
        'hbB_hbS_or_absolute_swirl_amplitude_materialized','first1_H_mean_or_finite_N_fields_zeroed')
    require(all(r['contract'][k] is True for k in true_keys) and all(r['contract'][k] is False for k in false_keys),'Construction contract changed')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    return {'passed':True,'signed_views':3,'view_terms':180,'unchanged_references_per_view':176,
        'four_signed_source_slots':list(LEADING_SLOTS),'new_source_zero_slots':[38,40],
        'macro_first_switch_leading_packets':3,'original_archives_target_reserved_six_and_N_retained':True,
        'global_terminal_closure':False}


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,expected in r[key].items():require(sha(base/name)==expected,'Changed pin: '+name);count+=1
    return count


def run():
    started=time.perf_counter();fixture=new_controls();o=OriginalFirstSwitchTerminalM15();r=o.build()
    code=dict(o.implementation_hashes);code.update({n:sha(ROOT/n) for n in CODE})
    r.update(input_hashes=o.input_hashes,implementation_hashes=code)
    r['checks']={'passed':True,'new_exact_zero_controls':fixture,'packet_invariants':check_packet(r,o.parent),'no_accepted_ancestor_replay':True}
    r['execution_seconds']=round(time.perf_counter()-started,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] and r['checks']['new_exact_zero_controls']['passed'],'Saved report failed')
    return {'passed':True,'pins':pins(r),'packet_invariants':check_packet(r,parent()),
        'serialized_only':True,'owner_or_fixture_or_ancestor_run':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-report',action='store_true');g.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
