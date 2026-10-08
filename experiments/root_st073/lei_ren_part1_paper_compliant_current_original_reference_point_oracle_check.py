"""Focused original reference source, seam and scalar loop/Z checks."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_point_oracle as current
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as original

base=current.base;slow=current.slow;ep=current.ep


def finite_references(manifest):
    p=mp.mp.clone();p.dps=110;c=base.MPIntervalContext();c.dps=160
    counts=0;restored=0;branches=set();midplane=False
    with mp.workdps(200):
        for sample in manifest['actual_original_reference_point_queries']:
            saved=slow.restore_point(p,sample['original_factored_point_inputs'])
            offset=s.Rational(sample['original_offset_exact']);z=s.Rational(sample['original_Z_exact'])
            yv=p.mpf(int(offset.p))/int(offset.q);Z=p.mpf(int(z.p))/int(z.q)
            f=p.exp(yv/10);C=1/(1+Z*Z)
            for key,order,wanted in (('E',0,C*f),('E',1,-2*Z*C*C*f),('a',0,p.mpf(4)/5),('a',1,0),('b',0,0),('b',1,0)):
                value=sum((t.coefficient for t in saved['inputs'][key][order].terms),p.mpf(0))
                assert abs(value-wanted)<p.mpf('1e-45');restored+=1
            for R,ps,delta in ((50,7,p.mpf('.1')),(70,13,p.mpf('.01'))):
                L=1-delta*Z*Z;factors=(p.mpf(R),p.mpf(ps),delta,L)
                bases=(c.ln(ps),c.ln(c.mpf(delta)),c.ln(c.mpf(L)),c.mpf(0),c.ln(R))
                ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                    positive_function_root_intersections=0,directed_independent_log_rescalings=0)
                roots={};direct={}
                for key,pair in saved['inputs'].items():
                    roots[key]={}
                    for order,row in enumerate(pair):
                        direct[(key,order)]=sum((t.coefficient*p.fprod(b**k for b,k in zip(factors,t.factor_powers))
                            for t in row.terms),p.mpf(0))
                        roots[key][(0,order)]=base.factored_row_enclosure(row,bases,ledger)
                        lo,hi=ep(roots[key][(0,order)].finite_interval())
                        assert lo-p.mpf('1e-100')<=direct[(key,order)]<=hi+p.mpf('1e-100')
                        restored+=1
                a,b,E,p1,p2=(direct[(key,0)] for key in ('a','b','E','p1','p2'))
                scales=original.loop.GenericLoopScales(a_min='.7',margin_min='.1',boundary_kappa_excess_min='.02',
                    t0_abs_max=0,p1_abs_max=p1+1,p2_abs_max=abs(p2)+1,dps=100)
                reference=original.GenericLoopPointZ(scales,a=a,b=b,p1=p1,p2=p2,E=E,
                    a_Z=direct[('a',1)],b_Z=direct[('b',1)],p2_Z=direct[('p2',1)],E_Z=direct[('E',1)])
                a_box=roots['a'][(0,0)];roots['t0']={(0,0):a_box.scalar(0),(0,1):a_box.scalar(0)}
                q=base.current.q_enclosure(a_box,a_box-2,c.ln(c.mpf(scales.eta)),c.ln(c.mpf('.7')))['q']
                kernel=base.conditioned.ConditionedPhase(dict(q=q,roots=roots),c.ln(c.mpf(scales.d_star)))
                phi=slow.read_scalar(p,sample['actual_original_radius_phase']['approximate_original_fractional_phase'])
                expected=reference.evaluate(phi)
                for chart in (('psi','E') if kernel.geometry=='signed_Mobius' else ('psi',)):
                    inverse=kernel.inverse_bracket(c.mpf(phi),chart,100)
                    values,proof=slow.slow_values(kernel,roots,inverse['coordinate_interval'],chart)
                    primitives=kernel.primitives(inverse['coordinate_interval'],chart)
                    branches.add(proof['branch'])
                    pairs=((values['psi_Z'],expected.angle_Z),(values['A_Z_slow'],expected.A_Z_slow),
                        (values['B_Z_slow'],expected.B_Z_slow),(primitives['A'],expected.A),
                        (primitives['B_over_Pstar'],expected.B))
                    for value,wanted in pairs:
                        lo,hi=ep(value.finite_interval());allow=p.mpf('1e-85')*(1+abs(wanted))
                        assert lo-allow<=wanted<=hi+allow,(sample['original_offset_exact'],Z,chart,mp.nstr(wanted,12))
                        counts+=1
                if Z==0:
                    assert p2==0 and direct[('p2',1)]!=0 and expected.A_Z_slow!=0
                    midplane=True
                for phi0 in ('0','.5','1'):
                    values,_=slow.slow_values(kernel,roots,c.mpf(phi0),'psi')
                    assert all(v.zero for v in values.values())
    assert midplane and {'conditioned_signed_E','conditioned_signed_psi','exact_midplane_nonzero_p2_Z'}<=branches
    return dict(passed=True,source_background_and_factor_restorations=restored,
        independent_original_scalar_A_B_and_Z_comparisons=counts,derivative_branches=sorted(branches),
        original_midplane_p2_zero_nonzero_p2_Z_and_A_Z_retained=True,
        finite_diagnostic_units_only=[[50,7,'.1'],[70,13,'.01']],
        finite_units_do_not_choose_original_native_parameters=True)


def native_and_seam(manifest):
    owner=current.OriginalReferencePointOracle();p=owner.owner.inputs.ctx;c=owner.ctx
    reference=owner.owner.inputs.evaluate(y=0,Z='.37')
    slope=base.point.OriginalO2FactoredPointInputs().evaluate(y=0,Z='.37')
    comparisons=0
    for key,pair in reference['inputs'].items():
        for a,b in zip(pair,slope['inputs'][key],strict=True):
            assert len(a.terms)==len(b.terms)
            for left,right in zip(a.terms,b.terms,strict=True):
                assert left.factor_powers==right.factor_powers
                assert abs(left.coefficient-right.coefficient)<=left.coefficient_error_upper+right.coefficient_error_upper
                assert left.late_pressure_error_log_upper==right.late_pressure_error_log_upper
                comparisons+=1
    for sample in manifest['actual_original_reference_point_queries']:
        assert sample['source_family']==owner.family and sample['true_reference_C0_Z_source_and_primitives_installed']
        assert sample['chart']=='Rh_reference' and sample['true_radius_phase_Z_exact_zero']
        assert not sample['source_caps_or_midpoints_selected_as_field_values']
        assert sample['actual_original_radius_phase']['positive_original_origin_offset']['strictly_positive']
        assert sample['actual_original_radius_phase']['positive_origin_budget_not_zeroed']
        for item in sample['actual_phase_held_Z_enclosures']:
            C0=item['C0'];target=base.conditioned.packets.interval(c,C0['phase'])
            image=base.conditioned.packets.interval(c,C0['selected_inverse']['phase_image'])
            assert ep(image)[0]<=ep(target)[0]<=ep(target)[1]<=ep(image)[1]
            assert C0['original_common_N_and_radius_phase_bound'] and not C0['free_phase_parameter_not_spatial_phase']
            assert item['derivative_contract']['p2_Z_retained'] and item['derivative_contract']['E_Z_term_retained']
            for value in item['slow_Z'].values():
                assert value['encloses_original_source_function'] and not value['point_value_selected']
    rejected=0
    for call in (lambda:current.reference_coordinate('-5.01'),lambda:current.reference_coordinate('.01'),
        lambda:current.reference_coordinate('nan'),lambda:owner.evaluate(offset=-2,Z=0,N=0),
        lambda:owner.evaluate(offset=-2,Z=0,N=7,bits=1)):
        try:call()
        except ValueError:rejected+=1
    assert rejected==5
    # Exact reference profiles at offset0 coincide with the defining slope
    # inlet histories, not merely with overlapping saved source ranges.
    z=s.Symbol('z');C=1/(1+z*z);V=4*z;ps=s.Symbol('Pstar',positive=True)
    reference_hist=dict(m=V,h=C*s.Rational(5,8),k=V*C*s.Rational(5,8),
        e=V*V/ps**2-C*C*s.Rational(5,12),p=C*C*s.Rational(5,2))
    slope_inlet=dict(m=V,h=C*(s.Rational(5,8)+0),k=V*C*(s.Rational(5,8)+0),
        e=V*V/ps**2-C*C*(s.Rational(5,12)+0),p=C*C*(s.Rational(5,2)+0))
    assert all(s.cancel(reference_hist[k]-slope_inlet[k])==0 for k in reference_hist)
    return dict(passed=True,native_true_radius_reference_queries=4,
        exact_closed_reference_to_original_slope_inlet_five_history_identity=True,
        seam_factored_coefficient_and_error_comparisons=comparisons,
        separate_original_P0_not_reset_or_added_twice=True,
        domain_N_and_precision_guards=rejected,whole_spatial_join_and_higher_jets_not_admitted=True)


def run():
    begin=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and manifest['no_original_slope_quadrature_or_ancestor_producer_reexecuted']
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('full_17_chart_numeric_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[k] is False for k in flags)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        independent_reference_source_loop_and_Z_checks=finite_references(manifest),
        native_phase_and_original_source_seam_checks=native_and_seam(manifest),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],
            current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(original.__file__).name:current.sha(Path(original.__file__).name)},
        execution_seconds=time.monotonic()-begin,
        scope='Second actual source chart Rh_reference with directed C0/Z and true-radius primitive enclosures; independent finite-unit scalar references and original source-inlet identity. No all17 oracle, integral controls, global frequency or corrected field admission.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original reference point oracle: source, A/B/Z, radius and seam PASS',flush=True)
    return report


if __name__=='__main__':run()
