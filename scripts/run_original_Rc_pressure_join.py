"""Produce the new Rc source congruence and quiet O3 signed primitive only."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/'experiments/root_st073'
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_finite_N_signed_source import load,sha,digest,require
from openai_ns_reconstruction.original_Rc_pressure_join import (
    PARENT,PARENT_SHA256,build,build_quiet_primitive,quiet_exact_controls)
from openai_ns_reconstruction.original_Rc_profile_join import build as profile_build
from openai_ns_reconstruction.original_Rc_endpoint_budget import rebuild_Rc_endpoint_budget
from openai_ns_reconstruction.original_Rc_endpoint_view import build_views,materialize_Rc178_view

NAME='lei_ren_part1_paper_compliant_current_original_Rc_pressure_join.json'
CODE=('src/openai_ns_reconstruction/original_Rc_profile_join.py',
    'src/openai_ns_reconstruction/original_Rc_pressure_join.py',
    'src/openai_ns_reconstruction/original_Rc_pressure_ports.py',
    'src/openai_ns_reconstruction/original_Rc_endpoint_budget.py',
    'src/openai_ns_reconstruction/original_Rc_endpoint_view.py',
    'scripts/run_original_Rc_pressure_join.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Changed accepted source parent')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed accepted parent required');return p


def pins(r):
    n=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():require(sha(base/name)==h,'Changed source pin: '+name);n+=1
    return n


def budget_check(budget):
    import mpmath as mp
    from mpmath.ctx_iv import MPIntervalContext
    from mpmath.libmp import mpf_cmp
    iv=MPIntervalContext();iv.dps=250;count=0
    exact=lambda pair:iv.mpf([mp.make_mpf(tuple(v)) for v in pair])
    for sample in budget['samples']:
        radius=sum((exact(v) for v in sample['radius_atoms']),iv.mpf(0))
        for component in sample['components']:
            reference=exact(component['endpoint_log_absolute_upper_exact_mpf_tuples'])+2*radius
            bound=exact(component['weighted_log_bound_exact_mpf_tuples'])
            require(mpf_cmp(bound._mpi_[0],reference._mpi_[0])<=0 and
                mpf_cmp(reference._mpi_[1],bound._mpi_[1])<=0,'200-digit directed budget failed 250-digit source reassembly')
            count+=1
    require(count==6,'Six independent C0/C1 bound rows required')
    return dict(passed=True,directed_lower_endpoint_rows=6,source_reassembly_precision=250,
        producer_precision=200,source_enclosures_only=True,original_function_or_quadrature_executed=False)


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,'Changed parent')
    require(r['source_family']==p['source_family'],'Different original family')
    require(r['quiet179_discrepancy_views']==p['quiet179_discrepancy_views'] and
        r['inherited_m15_packet']==p['inherited_m15_packet'],'Original archive or accepted views changed')
    profile,join,quiet=(r[k] for k in ('Rc_profile_join','Rc_pressure_join','quiet_O3_signed_primitive'))
    require(all(v['passed'] is True and v['source_executed'] is False for v in (profile,join,quiet)),
        'Source-only proof missing')
    require(profile['whole_O3_power_quiet'] is True and
        profile['input_source_prefix_sha256']==join['input_source_prefix_sha256'],'Different source profiles')
    require(join['required_profile_proof_digest']==digest(profile) and
        quiet['required_profile_proof_digest']==digest(profile),'Profile premise changed')
    require(len(join['background_mean_identities'])==4 and len(join['full_mean_identities'])==4 and
        len(join['actual_endpoint_maps'])==16,'Incomplete mean Z3 input congruence')
    require(all(v['passed'] is True for v in join['AST'].values()) and
        join['actual_physical_pressure_ports']['passed'] is True and
        all(join['actual_physical_pressure_ports']['checks'].values()) and
        join['same_source_Rc_boundary_pair_exact_zero'] is True and
        join['oppositely_signed_pair_C0_C1']==[0,0] and
        join['common_boundary_individually_zero_claim'] is False,'Wrong pressure endpoint contraction')
    require(len(quiet['actual_zero_mean_integrands'])==8 and
        all(v['passed'] is True for v in quiet['AST'].values()) and
        quiet['exact_source_quadrature_eliminated'] is True and
        quiet['all_quadratic_and_inlet_history_terms_retained'] is True and
        quiet['incoming_pressure_not_reset'] is True,'Incomplete exact quiet primitive')
    budget=r['Rc_remaining_lower_endpoint_budget']
    require(budget['passed'] is True and budget['contraction_applied'] is False and
        budget['no_aggregate_positive_bound_subtraction'] is True,'Endpoint budget was subtracted as a value')
    directed_budget=budget_check(budget)
    archive=load(DIRECTORY/r['inherited_m15_packet']['report'])
    require(len(r['Rc178_discrepancy_views'])==3,'Three bounded axial diagnostics required')
    for view,sample,old,bound_sample in zip(r['Rc178_discrepancy_views'],archive['samples'],p['quiet179_discrepancy_views'],budget['samples']):
        require(view['Z']==sample['Z']==old['Z'] and view['parent_quiet179_view_sha256']==digest(old),
            'Accepted quiet179 view changed')
        proof=view['source_proof']
        require(proof['Rc_pressure_join_sha256']==digest(join) and proof['profile_join_sha256']==digest(profile) and
            proof['remaining_lower_source_sha256']==digest(quiet) and
            proof['budget_sample_sha256']==digest(bound_sample),'Endpoint source theorem digest changed')
        old_other=[v for v in old['entries'] if v['parent_reduced_index'] not in (176,179)]
        new_other=[v for v in view['entries'] if v['parent_reduced_index']!=176]
        require(new_other==old_other,'Unrelated signed or error term changed')
        require(view['archived_aggregate_endpoint_caps']==[dict(index=i,term=sample['correlated_reduced_m15_C0_C1_terms'][i])
            for i in (176,179)],'Original positive endpoint cap archive changed')
        from mpmath.libmp import mpf_cmp
        from functools import cmp_to_key
        top=max((tuple(v['weighted_log_bound_exact_mpf_tuples'][1]) for v in bound_sample['components']),
            key=cmp_to_key(mpf_cmp))
        leaf=next(v for v in view['entries'] if v['parent_reduced_index']==176)['term']
        require(leaf['log_scale_terms'][0]['source_enclosure_exact_mpf_tuples']==[list(top),list(top)] and
            all(v['source_enclosure_exact_mpf_tuples']==[[1,1,0,1],[0,1,0,1]] for v in leaf['coefficient_C0_C1']),
            'Remaining lower enclosure lost its directed upper or symmetric error coefficient')
        expanded=materialize_Rc178_view(sample['correlated_reduced_m15_C0_C1_terms'],view['entries'],
            sample['independent_target'],source_proof=proof)
        require(digest(expanded)==view['expanded_Rc178_digest'] and view['view_terms']==178 and
            view['all_other_177_terms_unchanged'] is True,'Active endpoint contraction changed')
    for proof in (join,quiet):
        require(proof['physical_scale_materialized'] is False and
            proof['bounds_or_midpoints_registered_as_values'] is False,'Selected a scale, midpoint or bound value')
        for node in proof['exact_local_expression_nodes']:
            if node['operation']=='original_function_reference':
                require(node['report_sha256']==r['input_hashes'][node['report']],'Unpinned function reference')
    require(all(v is False for v in r['limits'].values()),'Global physical layer claimed')
    return dict(passed=True,ordinary_Z_orders=[0,1],actual_background_mean_identities=4,
        actual_full_mean_identities=4,actual_endpoint_substitutions=16,quiet_zero_mean_integrands=8,
        exact_same_source_Rc_boundary_pair=True,exact_O3_source_primitive=True,
        active_discrepancy_view_terms=178,source_endpoint_pair_contracted=True,
        directed_lower_endpoint_budget_check=directed_budget,
        original_signed_function_values_numerically_evaluated=False,
        accepted_ancestor_or_quadrature_replayed=False,global_terminal_closure=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    profile=profile_build(ROOT,p);join=build(ROOT,p,profile)
    quiet=build_quiet_primitive(ROOT,p,profile,join);controls=quiet_exact_controls()
    budget=rebuild_Rc_endpoint_budget(ROOT);views=build_views(ROOT,p,profile,join,quiet,budget)
    r=copy.deepcopy(p);r.update(stage='actual_Rc_pressure_source_join_and_exact_quiet_O3_primitive',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        Rc_profile_join=profile,Rc_pressure_join=join,quiet_O3_signed_primitive=quiet,
        Rc_remaining_lower_endpoint_budget=budget,Rc178_discrepancy_views=views,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{name:sha(ROOT/name) for name in CODE}))
    r['contract']=dict(actual_source_congruence_before_pressure_pairing=True,
        quiet_O3_original_source_quadrature_replaced_by_exact_primitive=True,
        ordinary_C1_and_incoming_pressure_retained=True,quiet179_views_unchanged=True,
        Rc_aggregate_endpoint_caps_archived_and_replaced_by_remaining_lower_bound=True,
        all_volume_caps_retained=True,independent_target_not_subtracted_again=True,
        old_viscous_endpoint_contractions_not_repeated=True,original_function_values_numerically_installed=False,
        accepted_ancestor_or_quadrature_replayed=False)
    r['checks']=dict(passed=True,new_quiet_exact_controls=controls,packet_invariants=packet_check(r,p),
        accepted_ancestors_replayed=False)
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True and
        r['checks']['new_quiet_exact_controls']['passed'] is True,'Saved passed new checks required')
    return dict(passed=True,pins=pins(r),packet_invariants=packet_check(r,parent()),serialized_only=True,
        owner_fixture_or_quadrature_replayed=False)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--write-report',action='store_true');group.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=r['execution_seconds'],bytes=path.stat().st_size,sha256=sha(path))))
    else:print(json.dumps(inspect(path),indent=2))
