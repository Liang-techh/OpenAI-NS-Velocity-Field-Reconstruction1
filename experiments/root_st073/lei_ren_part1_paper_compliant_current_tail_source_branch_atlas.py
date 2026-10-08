"""Original tail branch atlas with exact endpoint-focused source partitions.

Queries the genuine chart-generic phase/inverse/density interface. No global
flat extrapolation, unresolved-cell erasure or selected cap values are used.
Original axial phase1 and transition phase0 collars remain explicit.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals as oracle_source

tail=oracle_source.tail;rc=oracle_source.rc;common=oracle_source.common
HERE,PREFIX,sha=oracle_source.HERE,oracle_source.PREFIX,oracle_source.sha
N=oracle_source.N;encode=oracle_source.encode
NAME=PREFIX+'current_tail_source_branch_atlas.json'
RECEIPT=PREFIX+'current_tail_source_branch_atlas_check.json'
GATE='original_axial_transition_endpoint_focused_phase_density_source_atlas_executed'


def original_edges(chart):
    if chart=='O2_axial':
        return [Fraction(0)]+[1-Fraction(1,10**j) for j in range(2,21)]+[Fraction(1)]
    if chart=='O3_slope_mu':
        return [Fraction(0)]+[Fraction(1,10**j) for j in range(20,1,-1)]+[Fraction(1)]
    raise ValueError('Original axial or transition chart required')


def atlas(oracle,chart):
    edges=original_edges(chart);rows=[];attempts=[];pending=[(a,b,0) for a,b in zip(edges,edges[1:])]
    # Refine only actual refusals, with an explicit finite source-query budget.
    # This bounded atlas is progress toward whole-chart coverage, not a claim
    # that depth2 must cover every smooth microscopic cutoff seam.
    refinements=0
    while pending:
        left,right,depth=pending.pop(0)
        got=oracle.query(chart,left,right,N);record=got['record']
        record.update(source_partition_depth=depth,original_endpoint_focused_partition=True)
        if got['values'] is None and depth<2 and refinements<12:
            middle=(left+right)/2;attempts.append(record);refinements+=1
            pending[0:0]=[(left,middle,depth+1),(middle,right,depth+1)]
        else:rows.append(record)
        print('Original tail atlas source:',chart,str(left),str(right),record['status'],flush=True)
    covered=[r for r in rows if r['status']=='enclosed']
    unresolved=[r for r in rows if r['status']!='enclosed']
    summary=dict(chart=chart,exact_original_native_domain=['0','1'],
        original_seed_partition_edges=list(map(str,edges)),actual_source_queries=len(rows)+len(attempts),
        accepted_partition_cells=len(rows),enclosed_source_cells=len(covered),unresolved_source_cells=len(unresolved),
        refined_refused_parent_cells=len(attempts),maximum_subdivision_depth=2,maximum_refined_parents=12,
        actual_enclosed_phase_geometries={geometry:sum(
            cell['original_phase_first_jet_source'].get('geometry')==geometry
            for row in covered for cell in row['actual_original_phase_inverse_and_density_source']['actual_spatial_signed_density_Z_cells'])
            for geometry in ('flat','small_r_series','signed_Mobius')},
        exact_unresolved_native_windows=[r['exact_native_endpoints'] for r in unresolved],
        all_original_native_points_covered_by_enclosed_or_explicit_unresolved_cells=True,
        complete_chart_numerical_integral_admitted=False,
        enclosing_cumulative_source_cover_not_substituted_for_actual_inverse=True)
    return dict(summary=summary,ordered_original_source_partition=rows,original_refused_parent_source_queries=attempts)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();prior=json.loads((HERE/oracle_source.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,oracle_source,prior['source_family'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(tail.baseline.build_owner(bridge)))
        oracle=oracle_source.GenericTailPhaseIntegrals(owner);archives=[];summaries=[]
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            oracle.Z=Z;charts={chart:atlas(oracle,chart) for chart in ('O2_axial','O3_slope_mu')}
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),original_source_charts=charts,
                original_axial_ordinary_y_not_native_selector_derivative=True,
                actual_radius_phase_and_original_inverse_not_replaced_by_caps=True,
                full_source_partition_preserved_even_when_cutoff_unresolved=True,
                numerical_complete_tail_or_Rc_targets_admitted=False)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_tail_source_branch_atlas_'+tag+'.json.gz'
            if len(compressed)>=100*1024*1024:raise ValueError('Split large lossless atlas')
            (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
            archives.append(dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),
                uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()))
            summaries.append(dict(exact_Z_range=list(Z),charts={chart:value['summary'] for chart,value in charts.items()}))
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Source closures disagree: '+name)
            hashes[name]=digest
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            actual_original_endpoint_focused_source_atlas_archives=archives,actual_source_atlas_summaries=summaries,
            source_seed_evidence=seed,input_hashes=hashes,execution_seconds=time.monotonic()-began,
            full_original_axial_transition_function_phase_coverage_certified=False,
            numerical_complete_tail_integrals_or_Rc_targets_or_closure_admitted=False,
            **dict.fromkeys(common.current.FLAGS,False),
            scope='Genuine original O2_axial/O3_slope_mu exact endpoint-focused source partitions on both strict-sign tiles; actual phase/inverse/density C0/Z records or explicit refusals. No complete numerical tail or terminal targets/controls/global N/stress/recursion/full NS.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
