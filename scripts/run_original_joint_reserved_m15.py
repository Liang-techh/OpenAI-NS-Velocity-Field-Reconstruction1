"""Produce the joint full-support m15 consumer or inspect saved data only."""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_joint_reserved_m15 import (
    PARENT,PARENT_SHA256,FULL,FULL_SOURCE,OriginalJointReservedM15,JOINT_OPERATION,
    contract_reserved_target,materialize_joint_discrepancy_view,negate_component_sum,
    sha,digest,load,mp,MPIntervalContext,
)

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_joint_reserved_m15_numeric.json'
CODE=('src/openai_ns_reconstruction/original_reserved_full_support_source.py',
      'src/openai_ns_reconstruction/original_joint_reserved_m15.py',
      'scripts/run_original_joint_reserved_m15.py')


def require(value,message):
    if not value:raise AssertionError(message)


def rational_tuple(v):
    sign,man,exp,bits=v
    require(sign in (0,1) and man>=0 and bits>=0 and man.bit_length()==bits,'Finite exact mpf required')
    return F((-1 if sign else 1)*man)*(F(2)**exp)


def endpoints(row):
    raw=row.get('source_enclosure_exact_mpf_tuples')
    if raw is None:raw=[row[k+'_exact_mpf_tuple'] for k in ('lower','upper')]
    return [rational_tuple(v) for v in raw]


def check_joint_interval(coeff,components):
    for k,row in enumerate(coeff):
        wanted=(-sum(endpoints(v[k])[1] for v in components.values()),
                -sum(endpoints(v[k])[0] for v in components.values()))
        lo,hi=endpoints(row)
        require(lo<=wanted[0]<=wanted[1]<=hi,'Joint directed interval misses saved component sum')
        require(hi-lo<=(wanted[1]-wanted[0])*(1+F(1,10**60)), 'Unexplained joint interval inflation')
    return True


def controls():
    """Independent finite algebra and exact endpoint arithmetic only."""
    passed=[]
    def check(name,value):require(value,name);passed.append({'name':name,'passed':True})
    matrix=[[F(v) for v in row] for row in ((1,2,1),(0,1,3),(2,-1,1))]
    a=[row[:]+[F(int(i==j)) for j in range(3)] for i,row in enumerate(matrix)]
    for i in range(3):
        pivot=next(j for j in range(i,3) if a[j][i]);a[i],a[pivot]=a[pivot],a[i]
        scale=a[i][i];a[i]=[v/scale for v in a[i]]
        for j in range(3):
            if j!=i:
                scale=a[j][i];a[j]=[x-scale*y for x,y in zip(a[j],a[i])]
    inverse=[row[3:] for row in a]
    action=lambda row:[sum(row[j]*inverse[j][k] for j in range(3)) for k in range(3)]
    full_response=action(matrix[2]);check('full_weighted_row_inverse_action',full_response==[0,0,1])
    # Nonconstant angular inputs and a huge common pulse exercise both C1 rows.
    z=F(2,7);theta0=(1+z*z,2*z);theta1=(z**3,3*z*z)
    pre=(F(1,10**20)+z/10**25,F(1,10**25));pulse=(F(10**80)+z*z,2*z)
    full=tuple(x+y for x,y in zip(pre,pulse))
    result=tuple(-sum(r*q[k] for r,q in zip(full_response,(theta0,theta1,full)))+pulse[k] for k in (0,1))
    for k in (0,1):check('full_support_shared_pulse_cancelled_ordinary_'+str(k),result[k]==-pre[k])
    partial=action([F(1,3),F(1,2),F(2,5)])
    direct=tuple(-sum(r*q[k] for r,q in zip(partial,(theta0,theta1,full)))+pulse[k] for k in (0,1))
    reduced=tuple(-partial[0]*theta0[k]-partial[1]*theta1[k]-partial[2]*pre[k]+(1-partial[2])*pulse[k] for k in (0,1))
    check('partial_row_keeps_angular_and_residual_pulse',direct==reduced and direct!=tuple(-v for v in pre))
    check('partial_pressure_coefficient_is_not_one',partial[2]!=1)
    c=mp.mp.clone();c.dps=80;iv=MPIntervalContext();iv.dps=110
    encode=lambda x:{'source_enclosure_exact_mpf_tuples':[list(v) for v in iv.mpf(x)._mpi_]}
    parts={name:[encode(pair) for pair in values] for name,values in (
        ('inner',(['.1','.2'],['-.3','.4'])),('lower',(['-.5','.7'],['.1','.6'])),
        ('preRp',(['-.05','.02'],['-.1','.9'])))}
    q=negate_component_sum(c,iv,parts,source_bound=True)
    check('independent_rational_endpoint_inclusion',check_joint_interval(q,parts))
    for bad in (False,1):
        try:negate_component_sum(c,iv,parts,source_bound=bad)
        except ValueError:check('reject_unproved_response_'+str(bad),True)
        else:raise AssertionError('Unproved common source accepted')
    entries=[{'parent_reduced_index':i,'name':'reserved_n1_own_swirl_once' if i==182 else str(i),
        'parent_term_digest':'fixed'} for i in range(183) if i not in (38,39,40)]
    joint={'parent_reduced_index':182,'parent_term_digest':'fixed'}
    for bad in (False,1):
        try:contract_reserved_target(entries,joint,source_bound=bad)
        except ValueError:check('reject_unproved_target_contraction_'+str(bad),True)
        else:raise AssertionError('Unproved target contraction accepted')
    try:contract_reserved_target(entries[:-1],joint,source_bound=True)
    except ValueError:check('reject_missing_reserved_source',True)
    else:raise AssertionError('Missing reserved source accepted')
    check('negative_target_subtraction_sign',-full[0]-(-pulse[0])==-pre[0])
    return {'passed':True,'control_count':len(passed),'controls':passed,'accepted_fixture_or_quadrature_replayed':False}


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Pinned terminal source changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Accepted terminal source required');return p


def packet_check(r,p):
    require(r['source_family']==p['source_family'] and r['parent_report_sha256']==PARENT_SHA256,'Parent/family changed')
    proof=r['reserved_full_support_source'];cached=r['cached_component_source'];joint=r['joint_reserved_minus_target_source']
    require(proof['passed'] is True and all(v['passed'] is True for v in proof['AST'].values()),'Source full-support identity unproved')
    require(proof['correlation_proved'] is True and proof['graph_identity']['passed'] is True and
        proof['archive_provenance']['passed'] is True,'Reserved/target source lineage unproved')
    require(proof['source_full_support_identity']['passed'] is True and proof['response_minus_target']['passed'] is True,
        'Full-support response or target contraction not source-bound')
    require(proof['source_full_support_identity']['response_source']['native_compact_x']=='2' and
        proof['source_full_support_identity']['response_source']['exact_angular_response']==[0,0,1], 'Wrong original selected support row')
    require(proof['response_minus_target']['ordinary_Z_orders']==[0,1] and proof['response_minus_target']['before_boxing'] is True,
        'Ordinary source correlation missing')
    require(proof['response_minus_target']['shared_unit']['U_Z_exact_zero'] is True and
        all(proof['physical_scale_guards'].values()),'Physical pressure scale/sign source guard failed')
    require(proof['numeric_value_available'] is False and proof['numeric_full_closure'] is False and
        proof['global_closure_claimed'] is False,'Support proof promoted to full numerical closure')
    require(cached['passed'] is True and all(v['passed'] is True for v in cached['AST'].values()),'Cached component source unproved')
    require(cached['report']==FULL and cached['source']==FULL_SOURCE,'Wrong cached source')
    require(cached['report_sha256']==r['input_hashes'][FULL] and cached['source_sha256']==r['input_hashes'][FULL_SOURCE],'Cached component pin changed')
    require(cached['component_bounds_are_function_values'] is False and cached['component_bounds_or_midpoints_selected'] is False,'Bound selected as value')
    require(cached['source_errors_and_amplitude_derivatives_preserved'] is True and cached['unit_Z_derivative_exact_zero'] is True,'Source error/derivative lost')
    full=load(DIRECTORY/FULL)['samples']['whole_Z'];components=full['pressure_component_enclosures']['m15']
    selected={name:components[name]['C1_enclosures'] if name!='preRp' else components[name] for name in ('inner','lower','preRp')}
    require(cached['normalized_components']==selected and cached['normalized_components_digest']==digest(selected),'Cached signed component rows changed')
    require(cached['Z_domain']==full['Z'] and joint['enclosure_domain']==full['Z'],'Changed axial enclosure domain')
    require(check_joint_interval(joint['coefficient_C0_C1'],selected),'Joint enclosure invalid')
    require(joint['identity']=='reserved-(-pulse)=-full+pulse=-(inner+lower+preRp)' and joint['ordinary_Z_orders']==[0,1],'Wrong signed ordinary C1 recipe')
    require(joint['selected_compact_x']=='2' and joint['exact_response']==[0,0,1],'Partial response substituted')
    require(joint['source_C1_pairs']==p['joint_m15_pressure_source']['source_function_bindings'],'Wrong composed source handles')
    require(joint['common_normalization']==p['joint_m15_pressure_source']['normalized_expression']['shared_unit'],'Different physical units')
    require(joint['support_proof_digest']==digest(proof) and joint['cached_component_source_digest']==digest(cached),'Joint proof digest changed')
    require(joint['source_identity_before_boxing'] and joint['independent_target_consumed_once'] and joint['enclosure_not_function_value'],'Pulse source correlation lost')
    require(joint['source_centers_resolved'] is False and joint['physical_scale_materialized'] is False and joint['accepted_ancestor_or_quadrature_replayed'] is False,'Unsupported value or source replay')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'Original ledger witness changed')
    require(r['inherited_forward180_views_digest']==digest(p['signed_m15_query_views']),'Parent180 views changed')
    require(r['inherited_leading_packets_digest']==digest(p['macro_first_switch_leading_viscous_packets']),'Signed leading packet changed')
    ledger=load(DIRECTORY/r['inherited_m15_packet']['report'])
    require(len(r['joint_discrepancy_views'])==3,'Three source queries required')
    for view,before in zip(r['joint_discrepancy_views'],p['signed_m15_query_views']):
        require(view['Z']==before['Z'],'Axial query changed')
        sample=next(v for v in ledger['samples'] if v['Z']==view['Z'])
        original=sample['correlated_reduced_m15_C0_C1_terms'];raw=sample['raw_m15_C0_C1_terms']
        require(len(original)==183 and len(raw)==213 and digest(original)==view['parent_original_reduced183_digest'] and digest(raw)==view['raw213_digest'],'Original archives changed')
        require(view['entries'][:-1]==before['entries'][:-1] and view['parent_forward180_entries_digest']==digest(before['entries']),'Forward source category changed')
        require(view['view_terms']==180 and view['unchanged_parent_entries']==179 and view['unchanged_original_references']==175,'Wrong discrepancy scope')
        require(view['joint_reserved_minus_target_leaves']==1 and view['forward_six_inner_lower_cancelled'] is False,'Forward cancellation or duplicate reserved source')
        last=view['entries'][-1]
        require(last['operation']==JOINT_OPERATION and last['source_contraction_digest']==digest(joint),'Joint contraction changed')
        require(last['reserved_and_target_provenance']=={'reserved':before['entries'][-1],
            'independent_target':before['independent_target'],'source_identity_proof_digest':digest(proof)},'Reserved/target provenance changed')
        require(last['term']['coefficient_C0_C1']==joint['coefficient_C0_C1'] and last['term']['log_scale_terms']==before['independent_target']['log_scale_terms'],'Joint physical coefficient/unit changed')
        require(last['term']['actual_function_binding']==joint,'Joint function binding changed')
        expanded=materialize_joint_discrepancy_view(original,view['entries'],view['independent_target'])
        require(digest(expanded)==view['expanded_joint180_digest'] and view['original_forward180_expansion_digest']==before['expanded_signed180_digest'],'Expanded discrepancy changed')
        require(expanded[:6]==original[:6] and view['six_inner_lower_digest']==before['six_inner_lower_digest'],'Six signed inner/lower terms changed')
        require(view['independent_target']==before['independent_target'] and view['independent_target_digest']==before['independent_target_digest'],'Independent target definition changed')
        require(view['selected_N_evidence_digest']==before['selected_N_evidence_digest'] and view['largest_unchanged_caps']==before['largest_unchanged_caps'],'Finite N/caps changed')
        require(view['external_target_subtraction_required'] is False and view['independent_target_consumed_once'] is True,'Target is missing or subtracted twice')
        require(view['global_terminal_closure'] is False,'Global closure introduced')
    for key,group in (('input_hashes',r['input_hashes']),('implementation_hashes',r['implementation_hashes'])):
        require(all(group.get(name)==h for name,h in proof[key].items()),'Unpublished source binding pin')
    true=('typed_full_support_joint_input_consumed','joint_source_enclosure_computed','shared_pulse_cancelled_as_function_before_boxing',
        'reserved_and_target_consumed_once','other179_parent_entries_unchanged','raw213_reduced183_archives_unchanged',
        'forward_six_inner_lower_and_all_finite_N_caps_retained','exact_N_and_physical_units_retained')
    false=('accepted_ancestor_or_fixture_replayed','physical_width_or_amplitude_materialized',
        'caps_midpoints_or_overlap_selected_as_values','independent_target_subtracted_twice')
    require(all(r['contract'][k] is True for k in true) and all(r['contract'][k] is False for k in false),'Consumer contract changed')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    return {'passed':True,'joint_views':3,'view_terms':180,'unchanged_parent_entries':179,'unchanged_original_references':175,
        'joint_reserved_minus_target_leaves':1,'source_target_consumed_once':True,'full_support_joint_C1_enclosure':True,
        'original_archives_six_signed_inner_lower_N_and_all_finite_N_caps_retained':True,'global_terminal_closure':False}


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,expected in r[key].items():require(sha(base/name)==expected,'Changed source pin: '+name);count+=1
    return count


def run():
    start=time.perf_counter();fixture=controls();o=OriginalJointReservedM15();r=o.build()
    r.update(input_hashes=o.input_hashes,implementation_hashes=dict(o.implementation_hashes,**{n:sha(ROOT/n) for n in CODE}))
    r['checks']={'passed':True,'independent_new_consumer_controls':fixture,'packet_invariants':packet_check(r,o.parent),
        'no_accepted_ancestor_replay':True}
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True and r['checks']['independent_new_consumer_controls']['passed'] is True,'Passed saved report required')
    return {'passed':True,'pins':pins(r),'packet_invariants':packet_check(r,parent()),
        'serialized_only':True,'owner_fixture_or_ancestor_replayed':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-report',action='store_true');g.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
