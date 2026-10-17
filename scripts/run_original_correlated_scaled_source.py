"""Direct same-history scaled source bridge and canonical signed collar primitive."""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/'experiments/root_st073'
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_correlated_scaled_source import (
    PARENT,PARENT_SHA256,build_correlated_scaled_source,build_core_pressure_identity,load,sha,digest,require)
from openai_ns_reconstruction.original_scaled_signed_density import correlated_scaled_density,run_exact_controls
from openai_ns_reconstruction.original_signed_density_expression import SignedDensityExpressionGraph

NAME='lei_ren_part1_paper_compliant_current_original_correlated_scaled_source.json'
CODE=('src/openai_ns_reconstruction/original_scaled_signed_density.py',
    'src/openai_ns_reconstruction/original_correlated_scaled_source.py',
    'src/openai_ns_reconstruction/original_collar_signed_primitive.py',
    'src/openai_ns_reconstruction/original_compact_inlet_identity.py',
    'scripts/run_original_correlated_scaled_source.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted finite-N source parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    n=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():require(sha(base/name)==h,'Changed source pin: '+name);n+=1
    return n


def expression_controls(source):
    rows=[]
    for label,chart in source['charts'].items():
        graph=SignedDensityExpressionGraph();graph.nodes=chart['exact_local_expression_nodes'];report=chart['report']
        refs=[n for n in graph.nodes if n['operation']=='original_function_reference']
        bindings={(report,n['original_node']):F(n['original_node']%13+1,17) for n in refs}
        # Independent NEW rational source-port fixtures only. The number11 is
        # not an evaluation, approximation, or surrogate for the original N.
        bindings.update({(report,0):F(0),(report,1):F(1),(report,20):F(1,3),
            (report,25220):F(1,5),(report,18):F(11)})
        value=lambda n:bindings[(report,n)];ports=chart['original_signed_ports'];vy=ports['delta_V_y']
        dvy=([value(n) for n in vy['rows']] if vy['kind']=='exact_compact_delta_V_y_rows' else
            [value(a)/value(vy['N'])+value(b) for a,b in zip(vy['slow'],vy['fast'])])
        result=correlated_scaled_density([value(n) for n in ports['baseline_V']],
            [value(n) for n in ports['baseline_q']],[value(n) for n in ports['delta_V']],
            [value(n) for n in ports['delta_m']],dvy,F(1,3),F(1,5))
        values=graph.evaluate(bindings)
        for key,indices in chart['normalized_local_roots'].items():
            wanted=result[key];require([values[n] for n in indices]==[wanted.value,wanted.dz],
                'Direct same-history expression mismatch: '+label+'/'+key)
            rows.append(dict(name=label+'/'+key,passed=True))
    return dict(passed=True,controls=rows,control_count=len(rows),only_new_rational_fixtures=True,
        original_source_values_evaluated=False,fixture_N_is_not_original_N=True)


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,'Changed parent scope')
    require(r['quiet179_discrepancy_views']==p['quiet179_discrepancy_views'],'Quiet discrepancy views changed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'],'Original archives/provenance changed')
    require(r['source_family']==p['source_family'],'Construction family changed')
    source=r['correlated_scaled_source'];collar=r['canonical_collar_signed_primitive'];inlet=r['compact_inlet_source_identity']
    pressure=r['core_pressure_source_identity']
    require(source['passed'] is True and source['source_executed'] is False,'Direct recovery source unbound')
    require(all(v['passed'] is True for v in source['AST'].values()),'Changed direct recovery AST')
    for label,chart in source['charts'].items():
        require(len(chart['same_history_sum_proofs'])==7 and all(v['passed'] is True for v in chart['same_history_sum_proofs']),
            'Full and correction do not share exact histories')
        require(chart['full_minus_background_boxes_used'] is False and chart['caps_removed_or_narrowed'] is False,
            'Enclosure subtraction or cap replacement introduced')
        require(chart['nonzero_incoming_pressure_retained'] is True and chart['numerical_original_source_values_installed'] is False,
            'Original numerical values or pressure inlet changed')
        require(set(chart['normalized_local_roots'])=={'delta_q','delta_q_y','v1','v2','F1','F2'},'Scale sectors lost')
        nodes=chart['exact_local_expression_nodes']
        require(not any(n['operation']=='exp' for n in nodes),'Physical exponential materialized')
        for node in nodes:
            if node['operation']=='original_function_reference':require(node['report_sha256']==r['input_hashes'][node['report']],
                'Unpinned actual source reference')
    require(collar['passed'] is True,'Canonical signed collar primitive unbound')
    require(collar['source_executed'] is False,'Canonical constructor replayed')
    require(collar['delta_six_zero_claim'] is False,'Canonical singular sector promoted to six-source equality')
    require(collar['exact_controls']['passed'] is True and collar['caps_removed_or_narrowed'] is False,
        'Canonical singular controls or cap contract changed')
    require(inlet['passed'] is True and inlet['source_executed'] is False and inlet['owners_replayed'] is False,
        'Compact inlet same-source theorem unavailable')
    require(inlet['identity_count']==dict(correction=30,pressure_integrals=8,radius=1,total=39) and
        inlet['unresolved_dependencies']==[] and inlet['aggregate_finite_N_caps_retained'] is True,
        'Incomplete inlet normalization or cap removal')
    require(all(v['equal'] and v['same_actual_source_atom'] for v in
        inlet['correction']+inlet['pressure_integrals']+inlet['radius']),'Inlet not identical to original outer source')
    require(pressure['passed'] is True and pressure['true_forward_and_canonical_pressure_same_function'] is True and
        pressure['all_order_q_Taylor_identity']['passed'] is True,'Core pressure function identity unavailable')
    require(pressure['forward_point_approximations_equal_true_function_claim'] is False and
        pressure['flat_remainder_or_mass_enclosures_selected_as_values'] is False and pressure['delta_six_zero_claim'] is False,
        'Pressure approximation or function identity promoted to full six-source equality')
    require(r['contract']==dict(direct_signed_mean_history_used=True,Pstar_linear_quadratic_sectors_separate=True,
        incoming_pressure_and_all_finite_N_caps_retained=True,quiet179_views_unchanged=True,
        six_core_collar_equality_claimed=False,original_numerical_source_values_installed=False,
        accepted_ancestor_or_quadrature_replayed=False),'Changed stage contract')
    require(all(v is False for v in r['limits'].values()),'Global completion introduced')
    return dict(passed=True,actual_C5_charts=2,same_history_row_identities=14,Rc_inlet_function_identities=39,
        true_core_pressure_function_identity=True,all_order_q_Taylor_identity=True,ordinary_Z_orders=[0,1],
        original_source_values_evaluated=False,quiet179_views_and_all_caps_unchanged=True,global_terminal_closure=False)


def run():
    from openai_ns_reconstruction.original_collar_signed_primitive import build_signed_collar_primitive
    from openai_ns_reconstruction.original_compact_inlet_identity import build_compact_inlet_identity
    start=time.perf_counter();p=parent();pins(p)
    source=build_correlated_scaled_source(ROOT,p);collar=build_signed_collar_primitive(ROOT,p)
    inlet=build_compact_inlet_identity(ROOT,p)
    pressure=build_core_pressure_identity(ROOT,p)
    require(collar['passed'] is True,'Canonical singular source failed')
    algebra=run_exact_controls();require(algebra['passed'] is True,'Scaled density controls failed')
    expressions=expression_controls(source)
    r=copy.deepcopy(p);r.update(stage='direct_same_history_scaled_source_and_canonical_singular_primitive',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        correlated_scaled_source=source,canonical_collar_signed_primitive=collar,compact_inlet_source_identity=inlet,
        core_pressure_source_identity=pressure,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(direct_signed_mean_history_used=True,Pstar_linear_quadratic_sectors_separate=True,
        incoming_pressure_and_all_finite_N_caps_retained=True,quiet179_views_unchanged=True,
        six_core_collar_equality_claimed=False,original_numerical_source_values_installed=False,
        accepted_ancestor_or_quadrature_replayed=False)
    r['numerical_scope']=dict(normalized_arithmetic_and_directed_interval_backend_implemented=True,
        actual_original_signed_targets_and_interior_integrals_evaluated=False,
        current_source_graphs_are_functions_not_numerical_input_records=True,
        source_endpoint_equalities_are_function_identities_not_selected_box_values=True,
        physical_scale_and_exact_N_not_materialized=True)
    r['checks']=dict(passed=True,new_scaled_density_controls=algebra,new_expression_controls=expressions,
        packet_invariants=packet_check(r,p),accepted_ancestors_replayed=False)
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True and r['checks']['new_scaled_density_controls']['passed'] is True and
        r['checks']['new_expression_controls']['passed'] is True,'Saved passed checks required')
    return dict(passed=True,pins=pins(r),packet_invariants=packet_check(r,parent()),serialized_only=True,
        owner_fixture_or_quadrature_replayed=False)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-report',action='store_true');g.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=r['execution_seconds'],bytes=path.stat().st_size,sha256=sha(path))))
    else:print(json.dumps(inspect(path),indent=2))
