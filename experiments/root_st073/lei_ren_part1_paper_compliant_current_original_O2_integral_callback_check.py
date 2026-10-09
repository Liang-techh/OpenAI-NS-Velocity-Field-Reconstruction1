"""Fresh O2 graph dispatch, exact route/source cuts and own-rate checks."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_integral_callback as current

ep=current.ep;leaves=current.leaves


def interval(c,row):
    return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                  mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))


def overlaps(a,b):
    lo,hi=ep(a);left,right=ep(b);assert max(lo,left)<=min(hi,right),(lo,hi,left,right)


def native_frames(owner,manifest):
    dispatches=mass_checks=pieces=cuts=0;frames=[]
    for saved in manifest['actual_original_O2_integral_frames']:
        frame=owner.frame(Z=saved['original_Z_exact'],N=saved['explicit_candidate_N'],
            count=saved['accepted_defining_source_cell_level'],bits=saved['inverse_bits'])
        frames.append(frame);c=owner.c
        with mp.workdps(c.dps+40):
            previous={i:left for i,(left,right) in enumerate(owner.bounds)}
            totals={i:{key:c.mpf(0) for key in current.transport.RATES} for i in range(5)}
            direct={key:{jet:c.mpf(0) for jet in ('C0','Z')} for key in current.transport.RATES}
            for cell in frame.record['genuine_source_cell_phase_integral_records']:
                index=cell['route_index'];left,right=map(sy.Rational,cell['exact_y_cell'])
                a,b=owner.bounds[index];assert previous[index]==left and a<=left<right<=b;previous[index]=right
                source=cell['source'];sl,sr=map(sy.Rational,source['exact_y_cell'])
                assert sl<=left<right<=sr and sr-sl==sy.Rational(1,frame.count)
                assert source['source_family']==owner.family
                assert cell['accepted_whole_source_cell_restricted_as_function_range']
                assert cell['original_route_kernel_right_endpoint']==str(b)
                assert cell['original_coordinate_Jacobian_applied_once']==1
                assert cell['alternatives_unioned_not_added_as_source_masses']
                cuts+=int((left,right)!=(sl,sr))
                phase=cell['actual_original_phase_origin_and_increment']
                assert phase['explicit_candidate_N']==frame.N and phase['chart']=='O2_slope'
                assert phase['archived_N7_phase_rows_not_used'] and phase['exact_source_cell_phase_increment']=='N*(right-left)'
                for proof in cell['actual_source_piece_and_density_graph_records']:
                    assert proof['local_source_factors_canceled_before_finite_export']
                    assert proof['nonlinear_signed_DAG_before_source_and_phase_union']
                    graph=proof['actual_density_graph_root_contract']
                    assert graph['accepted_function_graph_executed_not_density_formula_rewritten']
                    assert graph['actual_C0_bound'] and graph['actual_phase_held_Z_bound']
                    assert proof['correlated_original_p2_Z_carrier']['every_original_coefficient_and_positive_late_pressure_error_retained']
                    pieces+=1
                for key,rate in current.transport.RATES.items():
                    mass=cell['exact_positive_own_rate_masses'][key];assert ep(mass)[0]>0;totals[index][key]+=mass
                    rr=c.mpf(rate.numerator)/rate.denominator
                    endpoint=c.mpf(int(b.p))/int(b.q);decay=c.exp(-rr*(1-endpoint))
                    for jet in ('C0','Z'):
                        direct[key][jet]+=cell['five_signed_C0_Z_density_hulls'][jet][key]*mass*decay
            assert all(previous[i]==b for i,(a,b) in enumerate(owner.bounds))
            p=mp.mp.clone();p.dps=c.dps+100
            for index,(a,b) in enumerate(owner.bounds):
                width=p.mpf(int((b-a).p))/int((b-a).q)
                for key,rate in current.transport.RATES.items():
                    r=p.mpf(rate.numerator)/rate.denominator
                    target=width if not rate else -p.expm1(-r*width)/r
                    lo,hi=ep(totals[index][key]);assert lo<=target<=hi;mass_checks+=1
            for role,row in owner.rows.items():
                contract=owner.require_integral(row);node=contract['integral_node']
                value=owner.dispatch(row,frame);assert value is frame.values[node]
                assert hasattr(value,'_mpi_') and not hasattr(value,'scale')
                original=interval(c,saved['actual_C0_Z_original_graph_integral_intervals'][str(node)])
                assert ep(value)==ep(original),(role,frame.count,frame.N)
                dispatches+=1
            for key in current.transport.RATES:
                for jet in ('C0','Z'):overlaps(direct[key][jet],frame.whole_values[key][jet])
    c=owner.c;point_checks=0
    with mp.workdps(c.dps+40):
        point=owner.provider.frame(chart='O2_slope',coordinate='.53',Z='.37',N=160)
        cell=next(row for row in frames[0].record['genuine_source_cell_phase_integral_records']
            if sy.Rational(row['exact_y_cell'][0])<sy.Rational(53,100)<sy.Rational(row['exact_y_cell'][1]))
        for key in current.transport.RATES:
            for jet in ('C0','Z'):
                value=owner.provider.source(owner.provider.rows['O2_slope','density_'+key+'_'+jet],coordinate='.53',Z='.37',N=160)
                overlaps(c.mpf(ep(value)),cell['five_signed_C0_Z_density_hulls'][jet][key]);point_checks+=1
        for node in owner.roles:
            overlaps(frames[0].values[node],frames[1].values[node])
    return dict(passed=True,fresh_original_O2_candidate_N_frames=len(frames),
        actual_issued_original_five_route_C0_Z_callbacks=dispatches,
        independent_positive_route_own_rate_total_mass_checks=mass_checks,
        certified_original_source_phase_and_derivative_pieces=pieces,
        exact_boundary_source_cell_restriction_cases=cuts,
        independent_current_point_source_DAG_and_whole_cell_enclosure_overlaps=point_checks,
        same_N_64_and_256_source_level_integral_enclosures_overlap=True,
        independent_cell_weighted_and_route_propagated_full_window_totals_overlap=True,
        archived_N7_phase_or_integral_rows_not_used=True,parent_producer_or_quadrature_suites_reexecuted=False)


def guards(owner):
    row=owner.rows[0,'m','C0'];source=owner.sources[owner.node_ids[0,'m','C0']]
    frame=owner.frame(Z='.37',N=160,count=64);g=owner.provider.built['graph'];before=leaves.digest_rows(g.nodes)
    contracts=[owner.require_integral(item) for item in owner.rows.values()]
    assert len(contracts)==50 and len({item['integral_node'] for item in contracts})==50
    assert len({tuple(item['exact_endpoints']) for item in contracts})==5
    rejected=0
    def reject(call):
        nonlocal rejected
        try:call()
        except (ValueError,TypeError,ArithmeticError,NotImplementedError):rejected+=1
        else:raise AssertionError('Invalid O2 integral request accepted')
    for call in (lambda:owner.require_integral(dict(row)),lambda:owner.require_integral(source),
        lambda:owner.dispatch(row,replace(frame)),lambda:owner.frame(Z=0,N=160),
        lambda:owner.frame(Z=2,N=160),lambda:owner.frame(Z='.37',N=159),
        lambda:owner.frame(Z='.37',N=True),lambda:owner.frame(Z='.37',N=160,count=100),
        lambda:owner.frame(Z='.37',N=160,bits=7),lambda:owner.parameter('logP'),
        lambda:owner.source(source,coordinate='.53',Z='.37',N=160)):
        reject(call)
    old=row['upper'];row['upper']=owner.route[-1]['upper']
    try:reject(lambda:owner.require_integral(row))
    finally:row['upper']=old
    old=source['source_node'];source['source_node']=old+1
    try:reject(lambda:owner.require_integral(row))
    finally:source['source_node']=old
    assert before==leaves.digest_rows(g.nodes)==owner.provider.graph_digest
    assert owner.mode not in ('original_source','synthetic_test')
    return dict(passed=True,exact_issued_density_root_kernel_endpoint_and_measure_bindings=50,
        invalid_row_frame_domain_frequency_grid_and_full_oracle_rejections=rejected,
        original_graph_unchanged=True,partial_fixed_nonzero_Z_provider_not_full_scalar_oracle=True)


def run():
    began=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('full_17_chart_source_or_24_cell_integral_oracle_installed',
        'full_scalar_control_evaluator_compatibility_installed','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2IntegralCallback()
    assert owner.family==manifest['source_family'] and owner.source_graph_sha256==manifest['source_graph_sha256']
    checks=dict(actual_O2_graph_integrals_and_exact_source_route_cuts=native_frames(owner,manifest),
        source_provenance_and_partial_scope_guards=guards(owner))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        source_graph_sha256=owner.source_graph_sha256,**checks,**dict.fromkeys(flags,False),
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Fresh original O2 five-route actual shared-N C0/Z integral callbacks, exact decimal endpoint/source-grid restrictions, original full density DAG roots, true phase and positive own-rate masses. Fixed nonzero-Z partial provider; no full-Z terminal proof, full evaluator/global N, actual recursion or corrected NS.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(result),indent=2).encode()+b'\n')
    print('Original O2 graph integral callbacks: exact cuts, source, phase and errors PASS',flush=True)
    return result


if __name__=='__main__':run()
