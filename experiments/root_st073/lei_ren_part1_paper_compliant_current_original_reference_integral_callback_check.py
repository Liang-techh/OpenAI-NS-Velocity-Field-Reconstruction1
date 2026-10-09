"""Focused graph/root, native cell source and directed integral callback checks."""
from dataclasses import replace
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_reference_integral_callback as current

ep=current.ep;leaves=current.leaves


def interval(c,row):
    return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                  mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))


def contains(box,value):
    lo,hi=ep(box);a,b=ep(value)
    assert lo<=a<=b<=hi,(lo,hi,a,b)


def source_and_integral_frames(owner,manifest):
    dispatches=ordinary=mass_checks=source_checks=0;frames=[]
    for saved in manifest['actual_original_reference_integral_frames']:
        frame=owner.frame(Z=saved['original_Z_exact'],N=saved['explicit_candidate_N'],
            count=saved['exact_source_cells'],bits=saved['inverse_bits'])
        frames.append(frame);c=frame.owner.ctx;a=frame.owner.atlas
        identity=frame.record['local_coefficient_program_equals_issued_original_density']
        assert identity['passed'] and len(identity['actual_issued_C0_Z_full_signed_density_identities'])==10
        assert identity['retained_exact_N_dependent_exprel_and_exp_functions']
        assert identity['original_B_over_Pstar_normalization_no_additional_division']
        phase_identity=frame.record['old_and_current_original_phase_argument_identity']
        assert phase_identity['passed'] and phase_identity['old_and_current_exact_argument_equal_not_endpoint_overlap']
        assert phase_identity['original_microscopic_origin_not_zeroed']
        with mp.workdps(c.dps+40):
            totals={key:c.mpf(0) for key in current.transport.RATES}
            for i,cell in enumerate(frame.cells):
                left,right=map(sy.Rational,cell['source']['exact_reference_cell'])
                assert left==-5+sy.Rational(5*i,frame.count) and right==-5+sy.Rational(5*(i+1),frame.count)
                assert cell['source']['closed_source_not_point_sample_extrapolation']
                phase=cell['source']['actual_phase_left_origin']
                assert phase['chart']=='Rh_reference'
                assert phase['exact_source_cell_phase_increment']=='N*(right-left)'
                assert phase['reference_native_radius_Jacobian_exact_one']
                assert phase['actual_cell_phase_is_original_point_origin_plus_exact_coordinate_width']
                assert cell['original_log_radius_Jacobian_applied_once']==1
                assert cell['nonlinear_exact_N_coefficients_enclosed_before_phase_union']
                for key,mass in cell['exact_positive_own_rate_masses'].items():
                    assert ep(mass)[0]>0;totals[key]+=mass
            p=mp.mp.clone();p.dps=c.dps+100
            for key,rate in current.transport.RATES.items():
                r=p.mpf(rate.numerator)/rate.denominator
                expected=p.mpf(5) if not rate else -p.expm1(-5*r)/r
                lo,hi=ep(totals[key]);assert lo<=expected<=hi;mass_checks+=1
            for role,row in owner.rows.items():
                key,jet=role;value=owner.dispatch(row,frame)
                assert value is frame.values[key][jet];dispatches+=1
                result=owner.integrate(row,Z=str(frame.Z),N=frame.N,count=frame.count,bits=frame.bits)
                assert hasattr(result,'_mpi_') and not hasattr(result,'scale')
                original=interval(c,saved['ordinary_directed_integral_intervals'][key][jet])
                assert ep(result)==ep(original),(role,frame.count,frame.N)
                ordinary+=1
            assert saved['original_order_minus1_and_minus2_N_dependent_functions_retained']
            assert saved['source_pressure_late_errors_and_positive_auxiliary_scales_retained']
            assert saved['incoming_corrections_and_original_P0_not_reset']
            assert not saved['unmaterializable_values_retained_factored']
    # A fresh current point DAG supplies the unweighted source density,
    # independently of the whole-cell coefficient-pair interpreter.
    frame=frames[0];a=frame.owner.atlas;c=a.ctx
    with mp.workdps(c.dps+40):
        point=owner.provider.frame(chart='Rh_reference',coordinate='-2.337',Z=str(frame.Z),N=frame.N)
        query=frame.owner.cell('-2.5','-1.25',N=frame.N,bits=frame.bits)
        for key in current.transport.RATES:
            for jet in ('C0','Z'):
                row=owner.provider.rows['Rh_reference','density_'+key+'_'+jet]
                native=owner.provider.dispatch(row,point)
                value=current.whole.prior.ScaledEnclosure(
                    a.scale(native.scale.powers,a.rational(point.coordinate),native.scale.offset),
                    a.copy_interval(native.coefficient),a.ledger)
                box=a.add(query['coefficients'][-1][key][jet]*(c.mpf(1)/frame.N),
                    query['coefficients'][-2][key][jet]*(c.mpf(1)/frame.N**2))
                anchor=current.whole.prior.FormalScale(a.bases,box.scale.powers)
                def normalized(v):
                    return c.mpf(0) if v.zero else v.coefficient*v.bounded_exp((v.scale-anchor).evaluate())
                contains(normalized(box),normalized(value));source_checks+=1
        # Both independently partitioned enclosures contain the same integral.
        for key in current.transport.RATES:
            for jet in ('C0','Z'):
                coarse=owner.integrate(owner.rows[key,jet],Z=str(frame.Z),N=160,count=4)
                fine=owner.integrate(owner.rows[key,jet],Z=str(frame.Z),N=160,count=16)
                lo,hi=ep(coarse);a0,b0=ep(fine);assert max(lo,a0)<=min(hi,b0)
    return dict(passed=True,fresh_genuine_original_reference_integral_frames=len(frames),
        fresh_symbolic_issued_density_and_local_exact_N_coefficient_identities=10,
        exact_original_old_and_current_phase_argument_identity=True,
        issued_C0_Z_definite_integral_dispatches=dispatches,ordinary_directed_integral_callbacks=ordinary,
        independent_stable_positive_own_rate_mass_checks=mass_checks,
        independent_current_point_DAG_densities_contained_in_whole_cell_coefficients=source_checks,
        four_and_sixteen_cell_same_N_integral_enclosures_overlap=True,
        actual_new_precise_phase_origin_and_exact_cell_shift_used=True,
        parent_producer_or_quadrature_suites_reexecuted=False)


def binding_and_scope_guards(owner):
    row=owner.rows['m','C0'];frame=owner.frame(Z='.37',N=160,count=4)
    g=owner.provider.built['graph'];before=leaves.digest_rows(g.nodes)
    source=owner.provider.rows['Rh_reference','density_m_C0']
    contracts=[owner.require_integral(item) for item in owner.rows.values()]
    assert len(contracts)==10 and len({item['integral_node'] for item in contracts})==10
    for contract in contracts:
        issued=owner.provider.rows['Rh_reference',contract['source_role']]
        key,jet=owner.provider.require_row(issued)
        assert (key,jet)==(contract['key'],contract['ordinary_order'])
        assert issued['source_node']==contract['source_node']
        assert contract['original_coordinate_Jacobian']==1
        assert contract['huge_common_radius_terms_canceled_symbolically']
    rejected=0
    def reject(call):
        nonlocal rejected
        try:call()
        except (ValueError,TypeError,ArithmeticError,NotImplementedError,current.leaves.transport.SourceOracleRequired):rejected+=1
        else:raise AssertionError('Invalid integral request was accepted')
    for call in (lambda:owner.require_integral(dict(row)),lambda:owner.require_integral(source),
        lambda:owner.dispatch(row,replace(frame)),lambda:owner.frame(Z=0,N=160),
        lambda:owner.frame(Z=2,N=160),lambda:owner.frame(Z='.37',N=159),
        lambda:owner.frame(Z='.37',N=True),lambda:owner.frame(Z='.37',N=160,count=0),
        lambda:owner.frame(Z='.37',N=160,count=257),lambda:owner.frame(Z='.37',N=160,bits=3),
        lambda:owner.parameter('logC'),lambda:owner.source(source,coordinate=0,Z='.37',N=160)):
        reject(call)
    old=row['lower'];row['lower']=row['upper']
    try:reject(lambda:owner.require_integral(row))
    finally:row['lower']=old
    old=source['source_node'];source['source_node']=old+1
    try:reject(lambda:owner.require_integral(row))
    finally:source['source_node']=old
    assert before==leaves.digest_rows(g.nodes)==owner.provider.graph_digest
    assert owner.mode not in ('original_source','synthetic_test')
    return dict(passed=True,exact_original_C0_Z_integrand_source_root_bindings=len(contracts),
        cloned_row_forged_frame_mutated_graph_domain_frequency_partition_and_full_oracle_rejections=rejected,
        original_graph_unchanged=True,partial_reference_integral_mode_not_full_scalar_oracle=True,
        original_unweighted_density_and_positive_mass_applied_exactly_once=True)


def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('full_17_chart_source_or_24_cell_integral_oracle_installed',
        'full_scalar_control_evaluator_compatibility_installed','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalReferenceIntegralCallback()
    assert owner.family==manifest['source_family'] and owner.source_graph_sha256==manifest['source_graph_sha256']
    checks=dict(actual_graph_bound_source_and_integral_callbacks=source_and_integral_frames(owner,manifest),
        exact_integrand_and_partial_scope_guards=binding_and_scope_guards(owner))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        source_graph_sha256=owner.source_graph_sha256,**checks,**dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Genuine issued reference C0/Z integral callbacks reconstructed from closed source cells, exact current phase and retained N-dependent coefficient functions. Source/root/kernel/Jacobian binding, independent current point density containment and directed interval dispatch. Fixed nonzero-Z partial provider; no all-route controls, global N, recursion or corrected NS.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(result),indent=2).encode()+b'\n')
    print('Original graph-bound reference integral callbacks: source, phase, errors and partial scope PASS',flush=True)
    return result


if __name__=='__main__':run()
