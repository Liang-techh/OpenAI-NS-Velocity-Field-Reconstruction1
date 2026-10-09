"""Scoped original tail-source, dlogx masses and affine connection checks."""
from dataclasses import replace
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_patch_terminal_integrals as current
import lei_ren_part1_paper_compliant_current_original_reference_integral_callback_check as reference_check

ep=current.ep;points=current.points;leaves=current.leaves


def coordinate(text):return s.E if text=='E' else points.terminal_coordinate(text)


def numeric(c,x):
    if x.is_Rational:return c.mpf(int(x.p))/int(x.q)
    q=points.log_coordinate(x);return c.exp(c.mpf(int(q.p))/int(q.q))


def log_gap(c,a,b):
    if a==b:return c.mpf(0)
    if a.is_Rational and b.is_Rational:
        q=(b-a)/a;return c.log1p(c.mpf(int(q.p))/int(q.q))
    return c.log(numeric(c,b))-c.log(numeric(c,a))


def source_and_masses(owner,manifest):
    frames=[];comparisons=mass_checks=microscopic=zero_cases=phase_checks=0
    for saved in manifest['actual_terminal_integral_frames']:
        a,b=map(coordinate,saved['exact_patch_x_interval'])
        frame=owner.subinterval(a,b,Z=saved['original_Z_exact'],N=saved['explicit_candidate_N'],count=saved['count'],bits=saved['bits']);frames.append(frame)
        c=frame.owner.ctx;p=mp.mp.clone();p.dps=c.dps+150
        knots=list(map(coordinate,saved['exact_partition_knots']))
        assert knots==owner.partition(a,b,saved['count']) and len(frame.cells)==len(knots)-1
        assert all(l<r for l,r in zip(knots,knots[1:]))
        assert frame.record['true_whole_cell_sources_not_point_quadrature']
        assert len(frame.record['local_coefficient_program_equals_issued_actual_patch_density']['actual_issued_C0_Z_full_signed_density_identities'])==10
        with mp.workdps(c.dps+40):
            totals={k:c.mpf(0) for k in leaves.transport.RATES}
            for i,item in enumerate(frame.cells):
                source=item['source'];l,r=map(coordinate,source['exact_patch_x_cell'])
                assert [l,r]==knots[i:i+2] and source['chart']=='actual_patch'
                assert source['closed_source_not_point_sample_extrapolation'] and source['nonlinear_original_coefficient_precedes_phase_union']
                phase=source['actual_phase_left_origin']
                assert phase['chart']=='actual_patch' and phase['explicit_candidate_N']==frame.N
                assert phase['source_graph_sha256']==owner.source_graph_sha256
                assert phase['exact_source_phase_increment']=='N*log(right/left)' and phase['original_patch_log_radius_measure']=='dx/x'
                lo,hi=ep(phase['actual_cell_phase_increment']);wanted=frame.N*log_gap(p,l,r)
                assert lo<=wanted<=hi;phase_checks+=1
                for k,rate in leaves.transport.RATES.items():
                    value=item['exact_positive_requested_endpoint_masses'][k];lo,hi=ep(value);assert lo>0
                    rr=p.mpf(rate.numerator)/rate.denominator;w=log_gap(p,l,r);d=log_gap(p,r,b)
                    wanted=w if not rate else p.exp(-rr*d)*(-p.expm1(-rr*w))/rr
                    assert lo<=wanted<=hi,(k,str(l),str(r));mass_checks+=1;totals[k]+=value
                    if r.is_Rational and l.is_Rational and r-l<s.Rational(1,10**900):microscopic+=1
            for k,rate in leaves.transport.RATES.items():
                rr=p.mpf(rate.numerator)/rate.denominator;w=log_gap(p,a,b)
                wanted=w if not rate else -p.expm1(-rr*w)/rr
                lo,hi=ep(totals[k]);assert lo<=wanted<=hi;mass_checks+=1
                lo,hi=ep(frame.decay[k]);assert lo<=p.exp(-rr*w)<=hi
                for j in ('C0','Z'):
                    original=reference_check.interval(c,saved['actual_ordinary_C0_Z_forcing'][k][j])
                    assert ep(original)==ep(frame.forcing[k,j]);comparisons+=1
            assert frame.record['unknown_incoming_correction_at_exact_left_endpoint_retained'] and frame.record['separate_original_P0_P0_Z_unchanged']
            if a==b:
                assert not frame.cells and all(ep(v)==(0,0) for v in frame.forcing.values())
                assert all(ep(v)==(1,1) for v in frame.decay.values());zero_cases+=1
    coarse,fine=frames[:2]
    for k in leaves.transport.RATES:
        for j in ('C0','Z'):
            a,b=ep(coarse.forcing[k,j]);l,r=ep(fine.forcing[k,j]);assert max(a,l)<=min(b,r)
    assert microscopic==5 and zero_cases==1
    return dict(passed=True,genuine_terminal_partial_integral_frames=len(frames),ordinary_C0_Z_integral_comparisons=comparisons,
        independent_higher_precision_positive_dlogx_mass_checks=mass_checks,
        independent_true_log_radius_phase_increment_checks=phase_checks,
        positive_microscopic_masses=microscopic,zero_width_identity_cases=zero_cases,
        exact_partition_no_gaps_or_duplicate_mass=True,coarse_fine_integral_overlap_consistency_only=True),frames


def independent_point_containment(owner,frame):
    a=frame.owner.atlas;c=a.ctx;p=owner.provider
    point=p.frame(coordinate='2.337',Z=str(frame.Z),N=frame.N)
    item=next(item for item in frame.cells if coordinate(item['source']['exact_patch_x_cell'][0])<=point.coordinate<=coordinate(item['source']['exact_patch_x_cell'][1]))
    l,r=map(coordinate,item['source']['exact_patch_x_cell']);query=frame.owner.cell(l,r,N=frame.N,bits=32);count=0
    with mp.workdps(c.dps+40):
        for k in leaves.transport.RATES:
            for j in ('C0','Z'):
                native=p.dispatch(p.rows['actual_patch','density_'+k+'_'+j],point)
                value=current.whole.prior.ScaledEnclosure(a.scale(native.scale.powers,points.patch_offset(c,point.coordinate),native.scale.offset),a.copy_interval(native.coefficient),a.ledger)
                box=a.add(query['coefficients'][-1][k][j]*(c.mpf(1)/frame.N),query['coefficients'][-2][k][j]*(c.mpf(1)/frame.N**2))
                anchor=current.whole.prior.FormalScale(a.bases,box.scale.powers)
                normalized=lambda v:c.mpf(0) if v.zero else v.coefficient*v.bounded_exp((v.scale-anchor).evaluate())
                reference_check.contains(normalized(box),normalized(value));count+=1
    return dict(passed=True,independent_true_point_signed_DAG_enclosures_contained_in_whole_source_coefficients=count,
        point_test_not_used_to_define_cell_source_or_integral=True)


def affine_connection(owner,history,frames,manifest):
    comparisons=0;counts=0
    for saved in manifest['actual_terminal_to_reference_O2_affine_extensions']:
        frame=next(f for f in frames if f.left==coordinate(saved['exact_patch_left']) and f.right==s.E and f.N==saved['explicit_candidate_N'])
        extended=owner.extend(frame,history,saved['reference_O2_right']);down=history.prefix(saved['reference_O2_right'],Z=str(frame.Z),N=frame.N);c=history.c
        with mp.workdps(c.dps+40):
            for k,rate in leaves.transport.RATES.items():
                rr=c.mpf(rate.numerator)/rate.denominator
                expected_decay=c.exp(-rr*((c.mpf(int(down.right.p))/int(down.right.q))+5+current.exact_log_width(c,frame.left,s.E)))
                a,b=ep(expected_decay);l,r=ep(extended['decay'][k]);assert max(a,l)<=min(b,r);comparisons+=1
                for j in ('C0','Z'):
                    expected=down.forcing[k,j]+down.decay[k]*history.current.ordinary(frame.forcing[k,j])
                    assert ep(expected)==ep(extended['forcing'][k,j])
                    original=reference_check.interval(c,saved['ordinary_C0_Z_forcing'][k][j]);assert ep(original)==ep(expected);comparisons+=1
        assert extended['record']['unknown_incoming_at_patch_left_retained'] and extended['record']['original_P0_P0_Z_not_reset']
        roles=extended['record']['actual_original_ten_Rh_incoming_root_roles'];assert len(roles)==10
        expected={str(pair[label]):dict(key=k,ordinary_order=j) for k,pair in owner.original_cell['outgoing'].items() for j,label in (('C0','value'),('Z','Z'))}
        assert roles==expected
        counts+=1
    H,I,J,D1,D2=s.symbols('H I J D1 D2')
    assert s.expand(D2*(D1*H+I)+J-(D2*D1*H+D2*I+J))==0
    return dict(passed=True,actual_terminal_to_reference_O2_history_compositions=counts,
        fresh_affine_forcing_and_multiplier_comparisons=comparisons,
        exact_original_ten_Rh_incoming_roots_are_actual_patch_outgoing=True,
        symbolic_affine_composition_identity=True,unknown_boundary_functions_not_selected_or_zeroed=True)


def guards(owner,history,frames):
    frame=frames[0];row=owner.provider.rows['actual_patch','density_m_C0'];g=owner.provider.built['graph'];before=leaves.digest_rows(g.nodes)
    calls=[lambda:owner.subinterval('1.7',2,Z='.37',N=160),lambda:owner.subinterval(2,'2.72',Z='.37',N=160),
        lambda:owner.subinterval('2.5',2,Z='.37',N=160),lambda:owner.tail(Z=0,N=160),lambda:owner.tail(Z=2,N=160),
        lambda:owner.tail(Z='.37',N=159),lambda:owner.tail(Z='.37',N=True),lambda:owner.tail(Z='.37',N=160,count=0),
        lambda:owner.tail(Z='.37',N=160,count=257),lambda:owner.tail(Z='.37',N=160,bits=3),
        lambda:owner.integrate(row,Z='.37',N=160),lambda:owner.extend(replace(frame),history,0),
        lambda:owner.extend(frames[3],history,0),lambda:owner.extend(frame,history,'1.01')]
    rejected=0
    for call in calls:
        try:call()
        except (ValueError,TypeError,ArithmeticError,NotImplementedError):rejected+=1
        else:raise AssertionError('Invalid tail/source/history request admitted')
    old=row['source_node'];row['source_node']=old+1
    try:
        try:owner.tail(Z='.37',N=160)
        except ValueError:rejected+=1
        else:raise AssertionError('Mutated original graph admitted')
    finally:row['source_node']=old
    assert leaves.digest_rows(g.nodes)==before
    return dict(passed=True,invalid_source_domain_Z_N_partition_precision_graph_frame_full_integral_and_composition_rejections=rejected,
        tail_integral_not_promoted_to_full_original_patch=True,exact_graph_unchanged=True)


@current.precise.phase.native.inlet.source_precision
def run():
    began=time.monotonic();manifest=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    history=current.partial.OriginalPartialHistoryFunctions()
    owner=current.OriginalPatchTerminalIntegrals(points.OriginalPatchTerminalPointLeaves(history.current.provider))
    checks,frames=source_and_masses(owner,manifest)
    connection=affine_connection(owner,history,frames,manifest)
    accepted=json.loads((current.HERE/current.partial.RECEIPT).read_bytes())
    assert accepted['all_passed'] and accepted[current.partial.GATE] and accepted['source_family']==owner.family
    for name,digest in accepted['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('full_original_patch_integral_installed','whole_Z_contract_installed','actual_source_owned_upstream_history_installed',
        'actual_five_controls_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(manifest[key] is False for key in flags)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        true_terminal_source_integrals_and_dlogx_masses=checks,
        independent_point_DAG_whole_source_containment=independent_point_containment(owner,frames[0]),
        actual_original_affine_tail_to_reference_O2_connection=connection,
        partial_tail_integral_and_history_guards=guards(owner,history,frames),**dict.fromkeys(flags,False),
        input_hashes={**owner.hashes,**accepted['input_hashes'],current.partial.RECEIPT:current.sha(current.partial.RECEIPT),
            current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(reference_check.__file__).name:current.sha(Path(reference_check.__file__).name)},
        execution_seconds=time.monotonic()-began,scope=manifest['scope'])
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(report),indent=2).encode()+b'\n')
    print('Original terminal patch: whole sources, dlogx integrals and affine history extension PASS',flush=True)
    return report


if __name__=='__main__':run()
