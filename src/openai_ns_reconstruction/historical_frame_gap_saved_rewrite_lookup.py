"""Resolve historical frame gap roles from the adapter's saved rewrite map.

Only direct key lookups into an authenticated, previously saved mapping are
used. No graph edges, source values, ancestor operators or new substitutions
are evaluated. Missing saved entries remain explicitly unresolved.
"""
from __future__ import annotations
import base64
import gzip
import hashlib
import json
from pathlib import Path

ARCHIVE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz'
ARCHIVE_SHA='442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
PRIOR='experiments/root_st073/lei_ren_part1_paper_compliant_current_historical_frame_gap_deferred_roles.json'
REPORT_SHA='caca47450da722fa290d5cb8d29d845e2c747b2d81b311317966ab38a2e33e25'
PRIOR_SHA='3ba2e27cb66459ca75592de04bf6359fa3295733cf87aa1b02f4efcf62d30ec6'


def _bytes(root: Path, name: str, digest: str) -> bytes:
    raw=(root/name).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('input hash mismatch: '+name)
    return raw


def recover(root: Path):
    report=json.loads(_bytes(root,PRIOR,REPORT_SHA))
    raw_result=base64.b64decode(report['artifacts']['tmp/historical_frame_gap_deferred_roles_v1/result.json']['base64'],validate=True)
    if hashlib.sha256(raw_result).hexdigest()!=PRIOR_SHA:raise ValueError('embedded prior result hash mismatch')
    prior=json.loads(raw_result)
    if prior.get('status')!='returned' or prior['result'].get('schema')!='historical_frame_gap_deferred_roles_v1':
        raise ValueError('prior direct-role packet is not the expected saved result')
    p=prior['result']
    raw=_bytes(root,ARCHIVE,ARCHIVE_SHA)
    archive=json.loads(gzip.decompress(raw))
    frame=archive['Rp_native_C5_pulse_source_bindings']
    if frame['input_substitution']!=p['archive_input_substitution']:
        raise ValueError('saved inlet map changed')
    mapping=frame['rewritten_nodes']
    if not isinstance(mapping,dict) or not mapping:raise ValueError('missing saved rewrite mapping')
    def lookup(node):
        if isinstance(node,bool) or not isinstance(node,int):raise ValueError('source row is not an opaque integer reference')
        key=str(node)
        if key not in mapping:
            return {'original_ref':node,'status':'not_in_saved_rewrite_map','rewritten_ref':None}
        value=mapping[key]
        if isinstance(value,bool) or not isinstance(value,int):raise ValueError('saved mapping target is not an opaque integer reference')
        return {'original_ref':node,'status':'present_in_saved_rewrite_map','rewritten_ref':value}
    def rows(nodes):
        if not isinstance(nodes,list) or len(nodes)!=4:
            raise ValueError('expected four ordinary-Z source rows')
        return [lookup(x) for x in nodes]
    old_global=p['original_global_m15_rows'];new_global=p['adapter_global_m15_rows']
    mapped_global=rows(old_global)
    if not isinstance(new_global,list) or len(new_global)!=4 or any(isinstance(v,bool) or not isinstance(v,int) for v in new_global) or any(q['rewritten_ref']!=v for q,v in zip(mapped_global,new_global)):
        raise ValueError('saved global pressure mapping does not match prior frame rows')
    e=p['original_gap_endpoint_roles']
    roles={'field_radius':lookup(p['original_own_gap_field_template']['radius']),
           'gap_m15':rows(p['original_gap_m15_rows']),
           'lower_R_q_y':rows(e['lower']['R_q_y']),
           'upper_R_q_y':rows(e['upper']['R_q_y']),
           'lower_q_plus_q_y':rows(e['lower']['q_plus_q_y']),
           'upper_q_plus_q_y':rows(e['upper']['q_plus_q_y'])}
    allrows=[roles['field_radius']]+[entry for k,v in roles.items() if k!='field_radius' for entry in v]
    return {'schema':'historical_frame_gap_saved_rewrite_lookup_v1',
      'status':'saved_rewrite_map_lookups_not_numerical_values',
      'prior_result_sha256':PRIOR_SHA,'frame_archive_raw_sha256':ARCHIVE_SHA,
      'map_owner':'Rp_native_C5_pulse_source_bindings.rewritten_nodes',
      'graph_namespace':'saved_Rp_native_C5_pulse_source',
      'global_original_to_frame_rows':mapped_global,
      'global_rows_match_saved_frame_pressure':True,
      'roles':roles,'roles_present_count':sum(x['status']=='present_in_saved_rewrite_map' for x in allrows),
      'roles_missing_count':sum(x['status']=='not_in_saved_rewrite_map' for x in allrows),
      'mapping_source':'saved adapter source map; direct lookups only',
      'new_substitution_executed':False,'source_graph_nodes_traversed':False,
      'old_constructors_called':False,'source_values_computed':False,
      'source_ids_are_numeric_values':False,'endpoint_tensor_complete':False,
      'local_to_global_m15_identity_proved':False,'P_equals_F_proved':False,'G_equals_F_proved':False,
      'actual_error_reduction_proved':False,'admissible_for_numeric_use':False,'full_goal_complete':False}
