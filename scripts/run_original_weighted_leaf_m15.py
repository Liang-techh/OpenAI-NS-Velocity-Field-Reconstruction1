"""Produce signed weighted leaves; inspect the saved overlay without replay."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_weighted_leaf_m15 import (
    PARENT,PARENT_SHA256,OriginalWeightedLeafM15,materialize_view,sha,digest,load,
)
DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_weighted_leaf_m15_numeric.json'
CODE=('src/openai_ns_reconstruction/original_weighted_leaf_kernel.py',
      'src/openai_ns_reconstruction/original_weighted_leaf_m15.py',
      'scripts/check_original_weighted_leaf_kernel.py','scripts/run_original_weighted_leaf_m15.py')


def require(value,message):
    if not value:raise AssertionError(message)


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Pinned weighted parent bytes changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'],'Passed weighted parent required');return p


def check_packet(r,p):
    require(r['parent_report_sha256']==PARENT_SHA256 and r['source_family']==p['source_family'],'Parent/family changed')
    require(r['source_guards']['passed'] and all(v['passed'] for v in r['source_guards']['AST'].values()),'Same-source leaf proof failed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'Ledger witness changed')
    require(r['inherited_weighted_queries_digest']==digest(p['weighted_remainder_queries']) and
            r['inherited_weighted_moment_digest']==digest(p['weighted_pulse_moments']),'Accepted queries/moments changed')
    require(len(r['signed_leaf_queries'])==19 and len(r['signed_m15_query_views'])==3,'Query/view count changed')
    regressions=0
    for row,old in zip(r['signed_leaf_queries'],p['weighted_remainder_queries']):
        require((row['Z'],row['phase'])==(old['Z'],old['phase']) and row['weighted_query_digest']==digest(old),'Weighted query order/source changed')
        pair=row['signed_pair'];require(pair['passed'] and pair['width_materialized'] is False and
            pair['source_rows_preserved_through_Z2'] and pair['same_source_identity_proved'],'Signed leaf pair failed')
        require([v['raw_index'] for v in pair['leaves']]==[46,47] and [v['width_power'] for v in pair['leaves']]==[2,3],'Original width powers lost')
        require(len(pair['operator_regressions'])==4 and all(v['passed'] and v['numerical_regression_not_function_identity'] for v in pair['operator_regressions']),'Operator regression failed')
        regressions+=len(pair['operator_regressions'])
        for leaf in pair['leaves']:
            require(len(leaf['source_ordinary_rows'])==3 and len(leaf['coefficient_C0_C1'])==2,'Missing ordinary Z2 or C1')
            require(all(v['approximation_is_exact_source_value'] is False and v['source_enclosure_is_not_a_selected_value'] for v in leaf['coefficient_C0_C1']),'Approximation promoted to truth')
    ledger_name=p['inherited_m15_packet']['report'];ledger=load(DIRECTORY/ledger_name)
    for view in r['signed_m15_query_views']:
        old=next(v for v in ledger['samples'] if v['Z']==view['Z'])
        raw=old['raw_m15_C0_C1_terms'];terms=old['correlated_reduced_m15_C0_C1_terms']
        require(view['raw213_digest']==digest(raw) and view['archived_reduced183_digest']==digest(terms),'Original archived ledger changed')
        require(view['view_term_count']==len(view['entries'])==183 and view['unchanged_parent_references']==181,'Wrong replacement scope')
        expanded=materialize_view(terms,view['entries'])
        require(digest(expanded)==view['expanded_signed_view_digest'],'Signed view expansion changed')
        replacements=[v for v in view['entries'] if v['operation']=='same_source_signed_leaf_replacement']
        require([v['raw_index'] for v in replacements]==[46,47] and [v['parent_reduced_index'] for v in replacements]==[41,42],'Wrong original source leaves replaced')
        query=next(v for v in r['signed_leaf_queries'] if v['Z']==view['Z'] and v['phase']=='1.0')
        require(view['query_pair_digest']==digest(query),'Wrong endpoint source query used in ledger')
        for entry in replacements:
            before=terms[entry['parent_reduced_index']];after=entry['term']
            leaf=next(v for v in query['signed_pair']['leaves'] if v['raw_index']==entry['raw_index'])
            require(after['coefficient_C0_C1']==leaf['coefficient_C0_C1'],'Ledger coefficient differs from signed source query')
            require(before['log_scale_terms']==after['log_scale_terms'],'Physical scale atom changed')
            binding=after['actual_function_binding'];require(binding['actual_signed_weighted_source_leaf'] and binding['error_enclosure_only'] is False,'Signed replacement still an unspecified cap')
            require(binding['same_source_identity_digest']==digest(r['source_guards']['identity']),'Wrong defining source identity')
        require(expanded[:6]==terms[:6] and digest(raw[:6])==view['six_inner_lower_digest'],'Six signed inner/lower changed')
        require(sum(v['name']=='reserved_n1_own_swirl_once' for v in expanded)==1,'Reserved contribution count changed')
        require(view['independent_target']==old['independent_target'] and view['independent_target_digest']==digest(old['independent_target']),'Independent target changed')
        require(view['independent_target']['actual_function_binding']['sign']==-1 and
                view['independent_target']['actual_function_binding']['physical_measure']=='-R*pressure_density*dR','P1 or a positive target substituted for independent negative m15')
        require(view['selected_N_evidence_digest']==digest(old['selected_N_evidence']),'Original common N changed')
        require(view['largest_unchanged_caps']==old['largest_remaining_individual_log_caps'] and view['closure_or_small_residual_claim'] is False,'Old dominant limitations changed')
        require(len(view['new_local_cap_comparisons'])==2,'Missing old/new coefficient budget comparison')
    require(r['contract']['accepted_ancestor_or_fixture_replayed'] is False and
            r['contract']['M0_W1_W2_or_core_recomputed'] is False and
            r['contract']['raw23_raw24_raw43_raw44_raw45_and_all_other_errors_retained'],'Unrelated work/ancestor changed')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    return {'passed':True,'signed_leaf_pairs':19,'operator_reconstructions':regressions,
        'Rp_signed_views':3,'view_terms':183,'unchanged_terms_per_view':181,'raw_leaves_replaced_in_views':[46,47],
        'old_raw213_reduced183_archives_unchanged':True,'independent_target_reserved_once_six_inner_lower_and_N_retained':True,
        'global_terminal_closure':False}


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,expected in r[key].items():require(sha(base/name)==expected,'Changed pin: '+name);count+=1
    return count


def run():
    started=time.perf_counter()
    from check_original_weighted_leaf_kernel import run_fixture
    fixture=run_fixture();owner=OriginalWeightedLeafM15();r=owner.build()
    code=dict(owner.implementation_hashes);code.update({name:sha(ROOT/name) for name in CODE})
    r.update(input_hashes=owner.input_hashes,implementation_hashes=code)
    r['checks']={'passed':True,'independent_new_leaf_fixture':fixture,
        'packet_invariants':check_packet(r,owner.parent),'no_accepted_ancestor_replay':True}
    r['execution_seconds']=round(time.perf_counter()-started,3)
    return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] and r['checks']['independent_new_leaf_fixture']['passed'],'Saved report failed')
    return {'passed':True,'pins':pins(r),'packet_invariants':check_packet(r,parent()),
            'serialized_only':True,'owner_or_fixture_or_ancestor_run':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-report',action='store_true');g.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
