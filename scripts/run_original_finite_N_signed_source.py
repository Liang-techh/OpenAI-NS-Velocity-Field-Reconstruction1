"""Build new signed C5 source expressions; inspect saved reports without replay."""
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
from openai_ns_reconstruction.original_finite_N_signed_source import (
    PARENT,PARENT_SHA256,build_signed_sources,load,sha,digest,require)
from openai_ns_reconstruction.original_inner_collar_source import build_inner_collar_source
from openai_ns_reconstruction.original_finite_N_density_kernel import (
    run_exact_controls,ordinary_jet,physical_fields,volume_density,endpoint_density)
from openai_ns_reconstruction.original_signed_density_expression import SignedDensityExpressionGraph

NAME='lei_ren_part1_paper_compliant_current_original_finite_N_signed_source.json'
CODE=('src/openai_ns_reconstruction/original_finite_N_density_kernel.py',
    'src/openai_ns_reconstruction/original_signed_density_expression.py',
    'src/openai_ns_reconstruction/original_finite_N_signed_source.py',
    'src/openai_ns_reconstruction/original_inner_collar_source.py',
    'scripts/run_original_finite_N_signed_source.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted quiet parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    n=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():require(sha(base/name)==h,'Changed source pin: '+name);n+=1
    return n


def expression_controls(packet):
    rows=[]
    for name,chart in packet['charts'].items():
        # A NEW exact rational port-shaped fixture, not values of these
        # original functions. Arbitrary jet rows test the expression adapter.
        graph=SignedDensityExpressionGraph();graph.nodes=chart['exact_local_expression_nodes']
        refs=[n for n in graph.nodes if n['operation']=='original_function_reference']
        report=chart['report'];bindings={(report,n['original_node']):F(n['original_node']%17+1,19) for n in refs}
        bindings.update({(report,0):F(0),(report,1):F(1),(report,20):F(1,3),
            (report,chart['delta']):F(1,5),(report,chart['Pstar']):F(2),
            (report,chart['native_variable']):F(3,2),(report,chart['physical_radius']):F(7)})
        value=lambda n:bindings[(report,n)]
        j=lambda source:ordinary_jet(*(value(n) for n in source[:3]))
        ports=chart['source_ports'];p=value(chart['Pstar'])
        b=physical_fields(j(ports['baseline_V'])*p,j(ports['baseline_normalized_Q'])*p)
        f=physical_fields(j(ports['full_physical']['U']),j(ports['full_physical']['Q']))
        z=ordinary_jet(F(1,3),F(1),F(0));de=F(1,5)
        direct=volume_density(f,z,de)-volume_density(b,z,de)
        direct_F=endpoint_density(f,z,de)-endpoint_density(b,z,de)
        L=1-z*z*de;D=1-z*z;Uy=j(ports['baseline_V_y'])*p
        Qy0=(z*(1+de)*b.U+2*z*Uy-D*b.Uz)/L-b.Q
        boundary=direct_F*49-(j(ports['full_physical']['Q_y'])-Qy0)*7
        roots=chart['signed_local_roots'];through=max(roots['delta_boundary_C0_C1'])
        values=graph.evaluate(bindings,through=through)
        for label,wanted in (('delta_v_C0_C1',direct),('delta_F_C0_C1',direct_F),('delta_boundary_C0_C1',boundary)):
            require([values[n] for n in roots[label]]==[wanted.value,wanted.dz],'Factored actual-port expression mismatch')
            rows.append(dict(name=name+'_'+label,passed=True))
    g=SignedDensityExpressionGraph();rejected=False
    try:g.reference('new_fixture','synthetic',0,source_nodes=[{'operation':'absolute_error_cap'}])
    except ValueError:rejected=True
    require(rejected,'Bound-as-value rejection failed');rows.append(dict(name='bound_as_original_function_rejected',passed=True))
    return dict(passed=True,controls=rows,control_count=len(rows),original_functions_evaluated=False,
        only_new_exact_rational_port_shaped_fixtures=True)


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,'Parent scope changed')
    require(r['source_family']==p['source_family'],'Different construction family')
    require(r['quiet179_discrepancy_views']==p['quiet179_discrepancy_views'],'Accepted discrepancy view changed')
    require(r['inherited_quiet179_digest']==digest(p['quiet179_discrepancy_views']),'Inherited view digest changed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'],'Original raw213/reduced183 provenance changed')
    core=r['inner_collar_source'];source=r['finite_N_signed_source']
    require(core['passed'] is True and core['gaps']==[] and core['source_executed'] is False,'Unbound inner/collar source')
    require(core['exact_identity']['delta_six_zero'] is False and core['numeric_value_available'] is False,'Unproved six-source equality introduced')
    for Z,rows in core['forward_six']['terms_by_Z'].items():
        require(len(rows['terms'])==6 and rows['all_term_errors_retained'] is True,'Six signed term errors lost')
        singular=rows['term_5_signed_source_remainder']
        require(singular['bound_is_not_a_value'] is True and singular['zero_claim'] is False,'Singular cap used as value')
    require(source['passed'] is True and all(v['passed'] is True for v in source['AST'].values()),'Density/FTC source changed')
    require(source['Rm_baseline_identical_graph_prefix']['passed'] is True,'Different Rm background prefix')
    require(set(source['charts'])=={'compact','Rm_patch'},'Actual source chart coverage changed')
    for label,chart in source['charts'].items():
        require(chart['common_exact_N']==18 and chart['extra_N_division'] is False,'Common N changed or divided twice')
        require(chart['incoming_pressure_difference_not_zeroed'] is True and chart['caps_removed_or_narrowed'] is False,'Pressure inlet or finite N error dropped')
        require(chart['numerical_signed_values_installed'] is False and chart['terminal_zero_claim'] is False,'Unevaluated source promoted to numerical result')
        require(len(chart['factored_six_volume_terms'])==6,'Quadratic density sectors lost')
        nodes=chart['exact_local_expression_nodes'];roots=chart['signed_local_roots']
        for key,value in roots.items():require(len(value)==2 and all(type(i) is int and 0<=i<len(nodes) for i in value),'Invalid C1 expression roots: '+key)
        for node in nodes:
            if node['operation']=='original_function_reference':
                require(node['report_sha256']==r['input_hashes'][node['report']],'Unpinned original source reference')
        for k,i in enumerate(roots['delta_volume_C0_C1']):
            require(nodes[i]['operation']=='definite_integral' and nodes[i]['physical_measure']=='R^2 dx/x' and
                nodes[i]['ordinary_Z_derivative_order']==k,'Wrong signed physical volume measure')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    require(r['contract']==dict(quiet179_views_and_all_caps_unchanged=True,six_signed_inner_lower_retained=True,
        target_and_quiet_already_consumed=True,actual_C5_function_ports_used=True,all_quadratic_terms_retained=True,
        full_F_and_viscous_endpoints_retained=True,incoming_pressure_not_reset=True,
        bounds_or_midpoints_selected_as_values=False,accepted_ancestors_replayed=False),'Changed source contract')
    return dict(passed=True,actual_source_charts=2,ordinary_Z_orders=[0,1],quiet179_views_unchanged=True,
        all_six_signed_inner_lower_and_finite_N_caps_retained=True,numerical_signed_values_installed=False,
        global_terminal_closure=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    source=build_signed_sources(ROOT,p)
    core=build_inner_collar_source(p,root=ROOT,directory=DIRECTORY)
    algebra=run_exact_controls();require(algebra['passed'] is True,'New exact algebra controls failed')
    # Fractions in controls are exact diagnostics, encoded without rounding.
    algebra={k:str(v) if isinstance(v,F) else v for k,v in algebra.items()}
    expressions=expression_controls(source)
    r=copy.deepcopy(p);r.update(stage='actual_finite_N_signed_density_and_inner_collar_sources',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        inherited_quiet179_digest=digest(p['quiet179_discrepancy_views']),
        finite_N_signed_source=source,inner_collar_source=core,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(quiet179_views_and_all_caps_unchanged=True,six_signed_inner_lower_retained=True,
        target_and_quiet_already_consumed=True,actual_C5_function_ports_used=True,all_quadratic_terms_retained=True,
        full_F_and_viscous_endpoints_retained=True,incoming_pressure_not_reset=True,
        bounds_or_midpoints_selected_as_values=False,accepted_ancestors_replayed=False)
    r['checks']=dict(passed=True,new_exact_density_controls=algebra,new_expression_controls=expressions,
        packet_invariants=packet_check(r,p),no_accepted_ancestor_replay=True)
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True,'Saved passed report required')
    require(r['checks']['new_exact_density_controls']['passed'] and r['checks']['new_expression_controls']['passed'],'Passed saved algebra receipts required')
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
