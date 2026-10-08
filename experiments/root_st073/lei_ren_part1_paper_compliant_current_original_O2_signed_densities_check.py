"""Independent changed-minus-original histories and Z differential checks."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_signed_densities as current
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as original

base=current.base;slow=current.slow;ep=current.ep


def finite_references(manifest):
    p=mp.mp.clone();p.dps=120;c=base.MPIntervalContext();c.dps=160
    es,vs=sy.symbols('E V');hist=current.recovery.history_densities(es,vs)
    jacobian={key:(sy.lambdify((es,vs),sy.diff(expr,es),'mpmath'),
        sy.lambdify((es,vs),sy.diff(expr,vs),'mpmath')) for key,expr in hist.items()}
    counts=0;comparisons=0;maximum=p.mpf(0)
    with mp.workdps(200):
        for sample in manifest['actual_original_O2_signed_density_queries']:
            saved=slow.restore_point(p,sample['actual_original_O2_factored_inputs'])
            Z=base.point.pressure.exact_Z(sample['original_Z_exact']);Z=p.mpf(int(Z.p))/int(Z.q)
            N=sample['explicit_candidate_N'];phi=slow.read_scalar(p,sample['actual_original_radius_phase']['approximate_original_fractional_phase'])
            for R,ps,delta in ((50,7,p.mpf('.1')),(70,13,p.mpf('.01'))):
                L=1-delta*Z*Z;factors=(p.mpf(R),p.mpf(ps),delta,L)
                bases=(c.ln(ps),c.ln(c.mpf(delta)),c.ln(c.mpf(L)),c.mpf(0),c.ln(R))
                ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                    positive_function_root_intersections=0,directed_independent_log_rescalings=0)
                roots={};direct={}
                for key,pair in saved['inputs'].items():
                    roots[key]={}
                    for order,row in enumerate(pair):
                        direct[(key,order)]=sum((term.coefficient*p.fprod(b**power for b,power in zip(factors,term.factor_powers))
                            for term in row.terms),p.mpf(0))
                        roots[key][(0,order)]=base.factored_row_enclosure(row,bases,ledger)
                a,b,E,p1,p2=(direct[(key,0)] for key in ('a','b','E','p1','p2'))
                EZ,V,VZ=direct[('E',1)],direct[('V',0)],direct[('V',1)]
                if sample['original_y_exact']=='1':roots['a'][(0,0)]=roots['a'][(0,0)].scalar(2)
                scales=original.loop.GenericLoopScales(a_min='.7',margin_min='.1',boundary_kappa_excess_min='.02',
                    t0_abs_max=0,p1_abs_max=p1+1,p2_abs_max=abs(p2)+1,dps=110)
                reference=original.GenericLoopPointZ(scales,a=a,b=b,p1=p1,p2=p2,E=E,
                    a_Z=direct[('a',1)],b_Z=direct[('b',1)],p2_Z=direct[('p2',1)],E_Z=EZ).evaluate(phi)
                a_box=roots['a'][(0,0)];roots['t0']={(0,0):a_box.scalar(0),(0,1):a_box.scalar(0)}
                q=base.current.q_enclosure(a_box,a_box-2,c.ln(c.mpf(scales.eta)),c.ln(c.mpf('.7')))['q']
                kernel=base.conditioned.ConditionedPhase(dict(q=q,roots=roots),c.ln(c.mpf(scales.d_star)))
                C0=kernel.evaluate(c.mpf(phi),bits=100);selected=C0['selected_inverse']
                primitives=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                jets,_=slow.slow_values(kernel,roots,selected['coordinate_interval'],selected['chart'])
                query=dict(kernel=kernel,roots=roots,ledger=ledger)
                # Use the same accepted graph from the source report owner;
                # graph dependencies are already hash checked before this call.
                graph=finite_references.graph
                got=current.BoundDensityGraph(graph,query,primitives,jets,N).outputs()
                factor=p.exp(reference.A/N);deltaE=E*p.expm1(reference.A/N);deltaV=reference.B/N
                EN=E*factor;VN=V+deltaV;ENZ=factor*(EZ+E*reference.A_Z_slow/N);VNZ=VZ+reference.B_Z_slow/N
                expected=dict(E_N=EN,V_N=VN,delta_E=deltaE,delta_V=deltaV,E_N_Z=ENZ,V_N_Z=VNZ,
                    delta_E_Z=ENZ-EZ,delta_V_Z=VNZ-VZ)
                before=current.recovery.history_densities(E,V);after=current.recovery.history_densities(EN,VN)
                differences={key:after[key]-before[key] for key in before}
                derivatives={key:fE(EN,VN)*ENZ+fV(EN,VN)*VNZ-fE(E,V)*EZ-fV(E,V)*VZ
                    for key,(fE,fV) in jacobian.items()}
                for actuals,target in ((got['velocities'],expected),(got['densities'],differences),(got['density_Z'],derivatives)):
                    for key,actual in target.items():
                        lo,hi=ep(actuals[key].finite_interval());allow=p.mpf('1e-90')*(1+abs(actual))
                        assert lo-allow<=actual<=hi+allow,(key,N,R,sample['original_Z_exact'])
                        maximum=max(maximum,lo-actual,actual-hi,p.mpf(0));comparisons+=1
                counts+=1
    return dict(passed=True,finite_original_changed_velocity_and_history_differential_cases=counts,
        independent_component_comparisons=comparisons,maximum_reference_outside_enclosure_discrepancy=maximum,
        reference_density_method='original history_densities(changed)-history_densities(original)',
        reference_Z_method='symbolic Jacobian of original history_densities with independent scalar loop-Z derivatives',
        finite_fixtures_do_not_select_native_parameters=True,no_defining_quadratures_or_ancestor_constructors_reexecuted=True)


def native_contracts(manifest):
    branches=set();tiny=False;midplane=False;cases=0
    for sample in manifest['actual_original_O2_signed_density_queries']:
        assert sample['source_family']==manifest['source_family']
        assert sample['cached_actual_defining_point_coefficients_and_errors_reused']
        assert sample['full_original_five_signed_density_values_and_Z_installed_on_this_O2_point']
        assert not sample['source_caps_or_midpoints_used_as_field_values']
        for row in sample['actual_signed_density_point_queries']:
            assert row['C0']['original_common_N_and_radius_phase_bound'] and not row['C0']['free_phase_parameter_not_spatial_phase']
            branches.add(row['phase_held_Z_contract']['branch']);contract=row['actual_graph_execution_contract']
            assert contract['accepted_function_graph_executed_not_density_formula_rewritten']
            assert contract['actual_C0_bound'] and contract['actual_phase_held_Z_bound']
            assert contract['tiny_original_factored_x_not_zeroed'] and contract['own_rates']==current.recovery.RATES
            assert contract['own_rate_units_with_common_S_Pstar']==current.recovery.UNITS
            assert set(row['five_signed_own_rate_density_enclosures'])==set(current.recovery.RATES)
            assert set(row['five_signed_own_rate_density_Z_enclosures'])==set(current.recovery.RATES)
            for values in (row['five_signed_own_rate_density_enclosures'],row['five_signed_own_rate_density_Z_enclosures']):
                assert all(value['encloses_original_source_function'] and not value['point_value_selected'] for value in values.values())
            if sample['original_y_exact']=='1':
                assert not row['modulated_normalized_velocity_and_Z_enclosures']['delta_E']['exact_zero']
                assert not row['five_signed_own_rate_density_enclosures']['p']['exact_zero'];tiny=True
            if sample['original_Z_exact']=='0':
                assert not row['five_signed_own_rate_density_Z_enclosures']['p']['exact_zero'];midplane=True
            cases+=1
        assert all(not sample[key] for key in ('actual_changed_five_moment_integral_evaluated',
            'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed','current_whole_N_selected'))
    assert cases==5 and tiny and midplane and len(branches)==4
    return dict(passed=True,actual_true_radius_signed_five_density_and_Z_queries=cases,
        tiny_endpoint_expm1_and_pressure_density_not_zeroed=True,nonzero_midplane_pressure_density_Z_preserved=True,
        accepted_graph_and_own_rate_units_preserved=True,no_point_samples_claimed_as_continuous_cell_integrals=True)


def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    c=base.MPIntervalContext();c.dps=260
    frame=base.FactoredGenericLoopScaleFrame(manifest['source_family'],c);finite_references.graph=frame.graph
    flags=('actual_changed_five_moment_integral_evaluated','numerical_original_source_point_or_integral_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        independent_scalar_changed_histories_and_Z_differentials=finite_references(manifest),
        actual_native_source_density_contracts=native_contracts(manifest),
        full_original_five_signed_density_values_and_Z_installed_on_actual_O2_points=True,
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began,
        scope='Original signed normalized O2 density and Z point enclosures, checked against independent changed-minus-original histories. Continuous-cell integration, own histories, global controls/N and recursion remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('True-radius original O2 signed densities: independent changed histories and Z differentials PASS',flush=True)
    return result


if __name__=='__main__':run()
