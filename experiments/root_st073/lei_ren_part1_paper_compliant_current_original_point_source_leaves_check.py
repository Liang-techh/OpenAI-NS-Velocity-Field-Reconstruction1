"""Focused point-root dispatch, native seams and independent density checks."""
from dataclasses import replace
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_point_source_leaves as current
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as original
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as functional

base=current.base;slow=current.slow;ep=current.ep


def finite_density_reference(owner,report):
    """New reference-chart densities checked with independent scalar loops."""
    p=mp.mp.clone();p.dps=120;c=base.MPIntervalContext();c.dps=160
    E,V=sy.symbols('E V');history=current.density.recovery.history_densities(E,V)
    jac={key:(sy.lambdify((E,V),sy.diff(expr,E),'mpmath'),sy.lambdify((E,V),sy.diff(expr,V),'mpmath'))
        for key,expr in history.items()}
    cases=0;comparisons=0;maximum=p.mpf(0)
    selected=[row for row in report['actual_original_point_source_frames']
        if row['chart']=='Rh_reference' or row['original_coordinate_exact']=='731/1000']
    with mp.workdps(200):
        for sample in selected:
            saved=slow.restore_point(p,sample['actual_defining_point_inputs_and_errors'])
            z=sy.Rational(sample['original_Z_exact']);z=p.mpf(int(z.p))/int(z.q)
            N=sample['explicit_candidate_N'];native_boxes=sample['actual_original_point_phase']['actual_phase_directed_boxes']
            # This exact literal probe checks the operator in finite units;
            # it is not substituted into the native source or its true phase.
            phi=p.mpf('.2317')
            for R,ps,delta in ((50,7,p.mpf('.1')),(70,13,p.mpf('.01'))):
                L=1-delta*z*z;factors=(p.mpf(R),p.mpf(ps),delta,L)
                bases=(c.ln(ps),c.ln(c.mpf(delta)),c.ln(c.mpf(L)),c.mpf(0),c.ln(R))
                ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                    positive_function_root_intersections=0,directed_independent_log_rescalings=0)
                roots={};direct={}
                for key,pair in saved['inputs'].items():
                    roots[key]={}
                    for order,row in enumerate(pair):
                        direct[key,order]=sum((term.coefficient*p.fprod(b**power for b,power in zip(factors,term.factor_powers))
                            for term in row.terms),p.mpf(0))
                        roots[key][0,order]=base.factored_row_enclosure(row,bases,ledger)
                a,b,e,p1,p2=(direct[key,0] for key in ('a','b','E','p1','p2'))
                ez,v,vz=direct['E',1],direct['V',0],direct['V',1]
                scales=original.loop.GenericLoopScales(a_min='.7',margin_min='.1',boundary_kappa_excess_min='.02',
                    t0_abs_max=0,p1_abs_max=p1+1,p2_abs_max=abs(p2)+1,dps=110)
                scalar=original.GenericLoopPointZ(scales,a=a,b=b,p1=p1,p2=p2,E=e,
                    a_Z=direct['a',1],b_Z=direct['b',1],p2_Z=direct['p2',1],E_Z=ez).evaluate(phi)
                aa=roots['a'][0,0];roots['t0']={(0,0):aa.scalar(0),(0,1):aa.scalar(0)}
                q=base.current.q_enclosure(aa,aa-2,c.ln(c.mpf(scales.eta)),c.ln(c.mpf('.7')))['q']
                kernel=base.conditioned.ConditionedPhase(dict(q=q,roots=roots),c.ln(c.mpf(scales.d_star)))
                inverse=kernel.evaluate(c.mpf(phi),bits=100)['selected_inverse']
                primitives=kernel.primitives(inverse['coordinate_interval'],inverse['chart'])
                jets,_=slow.slow_values(kernel,roots,inverse['coordinate_interval'],inverse['chart'])
                query=dict(kernel=kernel,roots=roots,ledger=ledger)
                values=current.density.BoundDensityGraph(owner.views[sample['chart']],query,primitives,jets,N).outputs()
                en=e*p.exp(scalar.A/N);vn=v+scalar.B/N
                enz=p.exp(scalar.A/N)*(ez+e*scalar.A_Z_slow/N);vnz=vz+scalar.B_Z_slow/N
                old=current.density.recovery.history_densities(e,v);new=current.density.recovery.history_densities(en,vn)
                target={key:new[key]-old[key] for key in old}
                targetZ={key:fE(en,vn)*enz+fV(en,vn)*vnz-fE(e,v)*ez-fV(e,v)*vz for key,(fE,fV) in jac.items()}
                for actuals,wanted in ((values['densities'],target),(values['density_Z'],targetZ)):
                    for key,value in wanted.items():
                        lo,hi=ep(actuals[key].finite_interval());allow=p.mpf('1e-85')*(1+abs(value))
                        assert lo-allow<=value<=hi+allow,(sample['chart'],sample['original_coordinate_exact'],key)
                        maximum=max(maximum,lo-value,value-hi,p.mpf(0));comparisons+=1
                cases+=1
    return dict(passed=True,new_reference_and_arbitrary_O2_scalar_density_cases=cases,
        independent_changed_minus_original_and_Z_Jacobian_comparisons=comparisons,
        maximum_reference_outside_enclosure_discrepancy=maximum,
        finite_diagnostic_units_only=[[50,7,'.1'],[70,13,'.01']],
        finite_probe_not_selected_as_native_source_phase_or_parameters=True)


def scalar_same_basis(value,bases,ledger):
    c=bases[0].ctx
    return current.prior.ScaledEnclosure(current.prior.FormalScale(bases,value.scale.powers,c.mpf(value.scale.offset)),
        c.mpf(value.coefficient),ledger)


def native_callbacks(owner,frames,report):
    dispatched=0;ordinary=0;formal=0;midplane=False;endpoint=False;fresh_O2=False
    with mp.workdps(350):
        for frame in frames:
            if frame.chart=='O2_slope' and frame.coordinate==sy.Rational(731,1000):
                # Fresh production evaluates defining functions here; the
                # checker may reuse their accepted coefficients plus errors.
                fresh_O2=True
            for key in current.transport.RATES:
                for order in ('C0','Z'):
                    row=owner.rows[frame.chart,'density_'+key+'_'+order]
                    point=owner.dispatch(row,frame)
                    assert point is frame.values[order][key] and point.ctx is frame.ctx
                    assert point.ledger is frame.query['ledger']
                    expected=owner.views[frame.chart]['five_signed_increment_rate_roots'][key] if order=='C0' else owner.views[frame.chart]['five_signed_increment_rate_first_derivatives']['Z'][key]
                    assert row['source_node']==expected
                    try:
                        numeric=owner.source(row,coordinate=str(frame.coordinate),Z=str(frame.Z),N=frame.N)
                        assert hasattr(numeric,'_mpi_') and not hasattr(numeric,'scale')
                        lo,hi=ep(numeric);assert lo<=hi and all(mp.isfinite(q) for q in (lo,hi));ordinary+=1
                    except ArithmeticError:
                        assert point is frame.values[order][key];formal+=1
                    dispatched+=1
            if frame.Z==0:
                assert frame.query['roots']['p2'][0,0].zero
                assert not frame.query['roots']['p2'][0,1].zero
                assert not frame.values['Z']['p'].zero;midplane=True
            if frame.chart=='O2_slope' and frame.coordinate==1:
                assert not frame.values['C0']['h'].zero and not frame.values['C0']['p'].zero;endpoint=True
            record=frame.record;phase=record['actual_original_point_phase']
            assert phase['source_phase_not_an_independent_free_angle'] and phase['original_width_Jacobians_unchanged']
            assert record['ordinary_Z_at_fixed_true_radius_phase'] and record['original_P0_P0_Z_not_modified']
            assert record['phase_union_is_one_point_error_cover_not_duplicate_source_mass']
            for piece in frame.pieces:
                assert piece['execution_contract']['B_over_Pstar_already_normalized_not_multiplied_again']
                assert piece['derivative']['p2_Z_retained'] and piece['derivative']['E_Z_term_retained']
    assert dispatched==report['original_role_bound_point_dispatches']==90
    assert ordinary==report['ordinary_interval_callbacks'] and formal==len(report['unresolved_ordinary_materializations_with_factored_source_retained'])
    assert midplane and endpoint and fresh_O2
    left=next(frame for frame in frames if frame.chart=='Rh_reference' and frame.coordinate==0)
    right=next(frame for frame in frames if frame.chart=='O2_slope' and frame.coordinate==0)
    assert left.Z==right.Z and left.N==right.N
    bases=left.query['kernel'].q.scale.bases;other=right.query['kernel'].q.scale.bases
    assert all(ep(a)==ep(b) for a,b in zip(bases,other,strict=True))
    seam=0
    for order in ('C0','Z'):
        for key in current.transport.RATES:
            a=left.values[order][key];b=scalar_same_basis(right.values[order][key],bases,a.ledger)
            difference=a-b
            assert ep(difference.coefficient)[0]<=0<=ep(difference.coefficient)[1],(order,key)
            seam+=1
    return dict(passed=True,actual_original_root_role_and_point_dispatch_checks=dispatched,
        actual_directed_ordinary_interval_callbacks=ordinary,unmaterializable_source_values_retained_factored=formal,
        original_reference_O2_common_point_C0_Z_seam_overlaps=seam,
        accepted_exact_reference_O2_source_inlet_identity_retained=True,
        nonzero_midplane_p2_Z_and_pressure_density_Z_retained=True,
        microscopic_original_endpoint_expm1_and_pressure_density_not_zeroed=True,
        new_arbitrary_O2_defining_point_with_errors_used=True)


def guards(owner,frames):
    first=frames[0];row=owner.rows[first.chart,'density_m_C0'];g=owner.built['graph'];before=current.digest_rows(g.nodes)
    calls=[lambda:owner.source_factored(dict(row),coordinate=str(first.coordinate),Z=str(first.Z),N=first.N),
        lambda:owner.dispatch(row,replace(first,N=first.N+1)),
        lambda:owner.dispatch(row,next(f for f in frames if f.chart!=first.chart)),
        lambda:owner.frame(chart='O3_power',coordinate=0,Z=0,N=257),
        lambda:owner.frame(chart='Rh_reference',coordinate=-6,Z=0,N=257),
        lambda:owner.frame(chart='O2_slope',coordinate='1.01',Z=0,N=257),
        lambda:owner.frame(chart='O2_slope',coordinate='.5',Z='1.1',N=257),
        lambda:owner.frame(chart='O2_slope',coordinate='.5',Z=0,N=159),
        lambda:owner.frame(chart='O2_slope',coordinate='.5',Z=0,N=True),
        lambda:owner.frame(chart='O2_slope',coordinate='.5',Z=0,N=257,bits=1),
        lambda:owner.source(row,coordinate=str(first.coordinate),Z=str(first.Z),N=first.N,phase='.3'),
        lambda:owner.parameter('Md'),lambda:owner.integrate(lambda x:x,0,1),
        lambda:functional.FunctionEvaluator(owner.built,oracle=owner,Z=0,N=257,ctx=first.ctx)]
    rejected=0
    for call in calls:
        try:call()
        except (ValueError,TypeError,NotImplementedError,current.transport.SourceOracleRequired):rejected+=1
    assert rejected==len(calls)
    old=row['source_node'];row['source_node']=old+1
    try:
        try:owner.require_row(row)
        except ValueError:rejected+=1
        else:raise AssertionError('Mutated original source graph accepted')
    finally:row['source_node']=old
    assert current.digest_rows(g.nodes)==before
    return dict(passed=True,invalid_root_frame_source_domain_Z_N_phase_and_full_oracle_rejections=rejected,
        exact_original_transport_graph_unchanged=True,partial_point_provider_not_promoted_to_full_scalar_oracle=True)


def run(owner=None,frames=None):
    began=time.monotonic();report=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    cached=owner is None
    if owner is None:
        owner=current.OriginalPointSourceLeaves();frames=[]
        for saved in report['actual_original_point_source_frames']:
            chart=saved['chart'];q=sy.Rational(saved['original_coordinate_exact']);z=sy.Rational(saved['original_Z_exact'])
            parent=owner.reference.owner if chart=='Rh_reference' else owner.O2.owner
            parent.inputs.point_cache[q,z]=slow.restore_point(parent.inputs.ctx,saved['actual_defining_point_inputs_and_errors'])
            frames.append(owner.frame(chart=chart,coordinate=str(q),Z=str(z),N=saved['explicit_candidate_N']))
    assert report[current.GATE] and report['source_family']==owner.family
    assert current.digest_rows(owner.built['graph'].nodes)==report['exact_original_transport_graph_sha256']
    flags=('full_17_chart_source_or_integral_oracle_installed','full_scalar_control_evaluator_compatibility_installed',
        'actual_original_numerical_integrals_evaluated','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(report[key] is False for key in flags)
    native=native_callbacks(owner,frames,report);rejections=guards(owner,frames);finite=finite_density_reference(owner,report)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        source_graph_sha256=owner.source_graph_sha256,actual_graph_bound_source_point_callbacks=native,
        independent_scalar_changed_histories_and_Z_density_derivatives=finite,point_scope_and_source_binding_guards=rejections,
        checked_defining_point_coefficients_and_error_cache_reused=cached,
        historical_producers_or_quadrature_suites_reexecuted=False,**dict.fromkeys(flags,False),
        input_hashes={**report['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(original.__file__).name:current.sha(Path(original.__file__).name),Path(functional.__file__).name:current.sha(Path(functional.__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Genuine Rh_reference/O2 graph-root C0/Z source point callbacks with defining errors, precise radius phase, native factors and directed ordinary tails. All-chart functions, parameter/integral callbacks, evaluated controls/global N and actual recursion remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(result),indent=2).encode()+b'\n')
    print('Original graph-bound reference/O2 point leaves: native dispatch, scalar density/Z and scope guards PASS',flush=True)
    return result


if __name__=='__main__':run()
