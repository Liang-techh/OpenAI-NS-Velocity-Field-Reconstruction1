"""Independent finite-unit restoration/loop references and native C0 scope.

The actual defining quadratures and ancestral computations are not rerun.
Diagnostic finite units test the new point/error basis attachment against
the original scalar loop. Native results retain their original factors.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_conditioned_primitives as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original


def read(c,record):
    if isinstance(record,dict) and 'exact_mpf_tuple' in record:return c.make_mpf(tuple(record['exact_mpf_tuple']))
    return c.mpf(record)


def row_from_record(c,record):
    terms=[]
    for term in record['terms']:
        late=term['late_pressure_error']['log_upper']
        terms.append(current.point.FactoredPointTerm(tuple(term['original_R_Pstar_delta_L_powers']),
            read(c,term['approximate_source_point_coefficient']),
            read(c,term['directed_finite_coefficient_absolute_error_upper']),None if late is None else read(c,late)))
    return current.point.FactoredPointRow(tuple(terms),record['original_y'],record['original_Z'])


def interval(c,record):return current.conditioned.packets.interval(c,record)


def finite_references(manifest):
    p=mp.mp.clone();p.dps=120;c=current.MPIntervalContext();c.dps=160
    restored=primitive_checks=0;discrepancies=[]
    with mp.workdps(200):
        for sample in manifest['actual_original_O2_conditioned_point_queries']:
            saved=sample['actual_original_O2_factored_inputs']
            Z=p.mpf(str(current.point.pressure.exact_Z(sample['original_Z_exact']).evalf(140)))
            rows={key:[row_from_record(p,row) for row in pair] for key,pair in saved['inputs'].items()}
            for R,ps,delta in ((50,7,p.mpf('.1')),(70,13,p.mpf('.01'))):
                L=1-delta*Z*Z;factors=(p.mpf(R),p.mpf(ps),delta,L)
                bases=(c.ln(ps),c.ln(c.mpf(delta)),c.ln(c.mpf(L)),c.mpf(0),c.ln(R))
                ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                    positive_function_root_intersections=0,directed_independent_log_rescalings=0)
                direct={};embedded={}
                for key,pair in rows.items():
                    for order,row in enumerate(pair):
                        # Independent direct restoration uses the original
                        # (R,Pstar,delta,L) order, not the formal basis order.
                        value=sum((term.coefficient*p.fprod(base**power for base,power in zip(factors,term.factor_powers))
                            for term in row.terms),p.mpf(0))
                        enclosure=current.factored_row_enclosure(row,bases,ledger)
                        lo,hi=current.point.endpoints(enclosure.finite_interval())
                        assert lo-p.mpf('1e-115')<=value<=hi+p.mpf('1e-115'),(key,order)
                        direct[(key,order)]=value;embedded[(key,order)]=enclosure;restored+=1
                if sample['original_y_exact']=='1':
                    embedded[('a',0)]=embedded[('a',0)].scalar(2)
                a,b,E,p1,p2=(direct[(key,0)] for key in ('a','b','E','p1','p2'))
                assert p1-2>p.mpf('.1') and a>p.mpf('.7') and b==0
                scales=original.GenericLoopScales(a_min='.7',margin_min='.1',boundary_kappa_excess_min='.02',
                    t0_abs_max=0,p1_abs_max=p1+1,p2_abs_max=abs(p2)+1,dps=110)
                reference=original.GenericShearLoop(scales,a=a,b=b,p1=p1,p2=p2,Utheta=E)
                a_box=embedded[('a',0)]
                q=current.current.q_enclosure(a_box,a_box-2,c.ln(c.mpf(scales.eta)),c.ln(c.mpf('.7')))['q']
                roots={key:{(0,0):embedded[(key,0)]} for key in ('a','b','E','p2')}
                roots['t0']={(0,0):a_box.scalar(0)}
                kernel=current.conditioned.ConditionedPhase(dict(q=q,roots=roots),c.ln(c.mpf(scales.d_star)))
                phi=read(p,sample['actual_original_radius_phase']['approximate_original_fractional_phase'])
                got=kernel.evaluate(c.mpf(phi),bits=80);assert got['status']=='enclosed'
                expected=reference.evaluate(phi)
                psi=current.point.endpoints(got['psi_fraction_interval'])
                expected_fraction=expected['angle']/(2*reference.ctx.pi)
                assert psi[0]-p.mpf('1e-95')<=expected_fraction<=psi[1]+p.mpf('1e-95')
                selected=got['selected_inverse'];values=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                for key,reference_key in (('A','A'),('B_over_Pstar','B')):
                    lo,hi=current.point.endpoints(values[key].finite_interval());actual=expected[reference_key]
                    assert lo-p.mpf('1e-95')<=actual<=hi+p.mpf('1e-95'),(key,R,sample['original_Z_exact'])
                    discrepancies.append(max(lo-actual,actual-hi,p.mpf(0)))
                primitive_checks+=1
    return dict(passed=True,independent_original_factor_unit_restoration_checks=restored,
        independent_original_scalar_loop_inverse_and_A_B_comparisons=primitive_checks,
        maximum_reference_outside_enclosure_discrepancy=max(discrepancies),
        finite_R_Pstar_delta_units_only=[[50,7,'.1'],[70,13,'.01']],
        scalar_reference_roundoff_allowance='1e-95',
        finite_units_and_fixture_scales_do_not_select_native_parameters=True,
        no_original_defining_quadratures_or_ancestor_constructors_reexecuted=True)


def native_contracts(manifest):
    counts={};maximum=mp.mpf(0);positive_tiny_endpoint=False
    c=current.MPIntervalContext();c.dps=260
    for sample in manifest['actual_original_O2_conditioned_point_queries']:
        assert sample['source_family']==manifest['source_family']
        assert sample['original_C0_inverse_and_primitives_enclose_true_O2_point_source']
        assert not any(sample[key] for key in ('scalar_A_B_point_values_selected','slow_Z_primitives_installed',
            'actual_changed_five_moment_integral_evaluated','numerical_original_source_point_or_integral_oracle_installed',
            'actual_five_controls_installed','current_whole_N_selected'))
        contract=sample['selected_scale_contract']
        assert contract['source_field_caps_not_selected_as_values'] and contract['selected_eta_dstar_not_materialized']
        assert contract['whole_original_input_margins_and_scale_choices_attached']
        basis=sample['source_factor_basis'];assert basis['order']==['logPstar','logdelta','logL','zero','logR']
        assert basis['same_original_R_Pstar_delta_L_bound'] and basis['original_positive_eta_delta_and_Pstar_inverse_sectors_retained']
        branch=sample['conditioned_source_geometry']['branch'];counts[branch]=counts.get(branch,0)+1
        for value in sample['actual_radius_phase_C0_inverse_and_A_B_enclosures']:
            assert value['status']=='enclosed' and value['inverse_installed_on_this_box']
            assert value['original_common_N_and_radius_phase_bound'] and not value['free_phase_parameter_not_spatial_phase']
            assert value['actual_O2_defining_point_inputs_and_errors_consumed']
            assert not value['source_caps_or_midpoints_used_as_field_values']
            image=interval(c,value['selected_inverse']['phase_image'])
            target=interval(c,value['phase']);il,ih=current.point.endpoints(image);tl,th=current.point.endpoints(target)
            assert il<=tl<=th<=ih
            width=ih-il;assert width<mp.mpf('1e-15');maximum=max(maximum,width)
            if sample['original_y_exact']=='1':
                assert basis['exact_a1_equals2_intersection_used']
                assert not value['geometry']['q']['exact_zero']
                assert all(not primitive['exact_zero'] for primitive in value['primitives'].values())
                positive_tiny_endpoint=True
    assert counts=={'signed_Mobius':2,'small_r_series':2} and positive_tiny_endpoint
    return dict(passed=True,true_original_source_radius_phase_C0_queries=4,source_geometry_counts=counts,
        maximum_directed_inverse_phase_image_width=maximum,
        positive_original_eta_q_and_endpoint_primitives_not_zeroed=True,
        exact_original_a1_equals2_avoids_rounded_excess=True,
        source_factor_and_error_contracts_preserved=True,all_global_installation_flags_false=True)


def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and manifest['native_source_parameter_recipes_checked_against_common_frame']
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    c=current.MPIntervalContext();c.dps=260
    scales=current.FactoredGenericLoopScaleFrame(manifest['source_family'],c)
    for key,record in manifest['selected_original_whole_input_scale_contract']['selected_positive_logs'].items():
        assert current.point.endpoints(interval(c,record))==current.point.endpoints(scales.logs[key])
    precisions=manifest['effective_numerical_precisions']
    assert precisions==dict(interval_digits=260,defining_point_coefficient_digits=50,
        directed_radial_coefficient_digits=90,radius_phase_digits=80,directed_radial_cells=4096,
        separate_precisions_with_all_errors_propagated=True)
    rejected=0
    for call in (lambda:current.OriginalO2ConditionedPrimitives(dps=100),
        lambda:current.OriginalO2ConditionedPrimitives(coefficient_dps=20),
        lambda:current.bits({'lower_exact_mpf_tuple':[0,1,0,1],'upper_exact_mpf_tuple':[0,2,0,2]})):
        try:call()
        except ValueError:rejected+=1
    assert rejected==3
    flags=('slow_Z_primitives_installed','actual_changed_five_moment_integral_evaluated',
        'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*current.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        independent_point_basis_and_scalar_loop_references=finite_references(manifest),
        actual_original_native_point_contracts=native_contracts(manifest),
        same_original_selected_scale_graph_bits_and_separate_precision_contract=True,
        precision_and_nonsingleton_selected_scale_guards=rejected,
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began,
        scope='Actual original source point/error and true-radius C0 inverse/A-B enclosures. Independent finite-unit scalar loop references test attachment without selecting native scales; slow-Z primitives, signed integral, whole oracle, global controls and recursion remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.encoded(result),indent=2)+'\n',encoding='utf8')
    print('True-radius O2 C0 primitives: independent source-factor/scalar-loop and native contracts PASS',flush=True)
    return result


if __name__=='__main__':run()
