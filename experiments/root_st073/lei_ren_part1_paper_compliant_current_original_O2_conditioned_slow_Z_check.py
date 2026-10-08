"""Focused independent original scalar-Z and actual native-source checks."""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_conditioned_slow_Z as current
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as original

base=current.base;ep=current.ep


def finite_references(manifest):
    p=mp.mp.clone();p.dps=120;c=base.MPIntervalContext();c.dps=160
    checks=charts=midplane=0;branches=set();maximum=p.mpf(0)
    with mp.workdps(200):
        for sample in manifest['actual_original_O2_slow_Z_point_queries']:
            saved=current.restore_point(p,sample['actual_original_O2_factored_inputs'])
            Z=base.point.pressure.exact_Z(sample['original_Z_exact'])
            Z=p.mpf(int(Z.p))/int(Z.q)
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
                if sample['original_y_exact']=='1':roots['a'][(0,0)]=roots['a'][(0,0)].scalar(2)
                scales=original.loop.GenericLoopScales(a_min='.7',margin_min='.1',boundary_kappa_excess_min='.02',
                    t0_abs_max=0,p1_abs_max=p1+1,p2_abs_max=abs(p2)+1,dps=110)
                reference=original.GenericLoopPointZ(scales,a=a,b=b,p1=p1,p2=p2,E=E,
                    a_Z=direct[('a',1)],b_Z=direct[('b',1)],p2_Z=direct[('p2',1)],E_Z=direct[('E',1)])
                a_box=roots['a'][(0,0)];roots['t0']={(0,0):a_box.scalar(0),(0,1):a_box.scalar(0)}
                q=base.current.q_enclosure(a_box,a_box-2,c.ln(c.mpf(scales.eta)),c.ln(c.mpf('.7')))['q']
                kernel=base.conditioned.ConditionedPhase(dict(q=q,roots=roots),c.ln(c.mpf(scales.d_star)))
                phi=current.read_scalar(p,sample['actual_original_radius_phase']['approximate_original_fractional_phase'])
                expected=reference.evaluate(phi)
                available=('psi','E') if kernel.geometry=='signed_Mobius' else ('psi',)
                for chart in available:
                    inverse=kernel.inverse_bracket(c.mpf(phi),chart,100)
                    values,proof=current.slow_values(kernel,roots,inverse['coordinate_interval'],chart)
                    branches.add(proof['branch']);charts+=1
                    for key,actual in (('psi_Z',expected.angle_Z),('A_Z_slow',expected.A_Z_slow),('B_Z_slow',expected.B_Z_slow)):
                        lo,hi=ep(values[key].finite_interval());allow=p.mpf('1e-95')*(1+abs(actual))
                        assert lo-allow<=actual<=hi+allow,(key,chart,R,sample['original_Z_exact'],mp.nstr(actual,15),mp.nstr(lo,15),mp.nstr(hi,15))
                        maximum=max(maximum,lo-actual,actual-hi,p.mpf(0));checks+=1
                if Z==0:
                    assert p2==0 and direct[('p2',1)]!=0 and expected.A_Z_slow!=0
                    midplane+=1
                for phase in ('0','.5','1'):
                    values,_=current.slow_values(kernel,roots,c.mpf(phase),'psi')
                    assert all(value.zero for value in values.values())
    assert midplane==2 and {'conditioned_signed_psi','conditioned_signed_E','exact_midplane_nonzero_p2_Z','regular_small_r_Fourier'}<=branches
    return dict(passed=True,independent_original_scalar_Z_component_comparisons=checks,
        inverse_coordinate_chart_comparisons=charts,derivative_branches=sorted(branches),
        exact_midplane_nonzero_p2_Z_references=midplane,maximum_reference_outside_enclosure_discrepancy=maximum,
        exact_periodic_and_half_period_Z_symmetry_checked=True,
        finite_diagnostic_units_only=[[50,7,'.1'],[70,13,'.01']],
        native_parameters_not_selected_by_finite_fixtures=True,
        no_defining_source_quadratures_or_ancestor_constructors_reexecuted=True)


def native_contracts(manifest):
    counts={};midplane=False;endpoint=False
    for row in manifest['actual_original_O2_slow_Z_point_queries']:
        assert row['source_family']==manifest['source_family'] and row['true_radius_phase_Z_exact_zero']
        assert row['slow_Z_primitives_installed_on_this_O2_point'] and row['cached_actual_defining_quadrature_bits_and_errors_reused']
        assert not row['source_caps_or_midpoints_used_as_field_values']
        assert all(not row[key] for key in ('actual_changed_five_moment_integral_evaluated',
            'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed','current_whole_N_selected'))
        assert row['source_factor_basis']['same_original_R_Pstar_delta_L_bound']
        for item in row['actual_phase_held_Z_enclosures']:
            C0=item['C0'];proof=item['derivative_contract'];counts[proof['branch']]=counts.get(proof['branch'],0)+1
            assert C0['original_common_N_and_radius_phase_bound'] and not C0['free_phase_parameter_not_spatial_phase']
            assert proof['p2_Z_retained'] and proof['E_Z_term_retained']
            assert proof['positive_implicit_denominator_lower']==1
            assert proof['positive_rho_s_hinv_not_rounded_to_zero']
            assert all(value['encloses_original_source_function'] and not value['point_value_selected'] for value in item['slow_Z'].values())
            if row['original_Z_exact']=='0':
                assert proof['branch']=='exact_midplane_nonzero_p2_Z' and not proof['u_Z']['exact_zero']
                assert not item['slow_Z']['A_Z_slow']['exact_zero'];midplane=True
            if row['original_y_exact']=='1':
                assert not row['source_geometry']['q']['exact_zero']
                assert not item['slow_Z']['A_Z_slow']['exact_zero'];endpoint=True
    assert midplane and endpoint and counts.get('conditioned_signed_E')==1
    return dict(passed=True,actual_original_radius_phase_slow_Z_queries=len(manifest['actual_original_O2_slow_Z_point_queries']),
        derivative_branch_counts=counts,midplane_p2_zero_with_nonzero_p2_Z_and_A_Z_preserved=True,
        positive_tiny_endpoint_q_and_Z_rows_retained=True,original_E_Z_product_term_retained=True,
        ordinary_point_Z_rows_consumed_once_no_extra_L_derivative=True,
        exact_global_frequency_and_full_integral_flags_remain_false=True)


def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and manifest['no_original_defining_quadratures_or_ancestor_constructors_reexecuted']
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    assert current.exact_derivative_identities()['passed'] and manifest['exact_original_derivative_identities']['passed']
    flags=('actual_changed_five_moment_integral_evaluated','numerical_original_source_point_or_integral_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        independent_original_scalar_Z_references=finite_references(manifest),actual_native_source_Z_contracts=native_contracts(manifest),
        exact_transformed_and_Fourier_derivative_identities=True,
        slow_Z_primitives_installed_on_actual_O2_points=True,**dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Local O2 source/error and true-radius phase-held derivative enclosures. Original scalar finite-unit references exercise both conditioned charts; native scale selection, density integrals, all-chart oracle, controls and recursion remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('True-radius O2 slow-Z: original scalar-Z and native factor contracts PASS',flush=True)
    return result


if __name__=='__main__':run()
