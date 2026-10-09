"""Check true signed inverse/Z evaluations, actual phase and local C1 integrals."""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions as current
import lei_ren_part1_paper_compliant_current_original_Rm_conditioned_phase_density_check as scalar_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):return current.encode(current.serialized(value))


def independent_signed_inverse():
    def tested(source,qr,dstar,phi):
        _,branches,_=current.signed.signed_u_branches(source,dstar)
        got=[current.branch_Z_functions(source,qr,dstar,phi,branch) for branch in branches]
        assert all(row['values'] is not None for row in got)
        return dict(values={name:current.union([row['values'][name] for row in got]) for name in current.OUTPUTS},
            record=dict(status='enclosed',geometry=got[0]['record']['geometry']))
    proxy=SimpleNamespace(**{**vars(scalar_check.current),'Z_FIRST':tested})
    with patch.object(scalar_check,'current',proxy),mp.workdps(170):
        cases=scalar_check.scalar_fixtures()
    return dict(passed=True,independent_original_scalar_inverse_Z_phase_and_exact_nonlinear_density_diagnostics=cases,
        actual_current_conditional_branch_Z_backend_used=True,
        finite_difference_and_scalar_references_diagnostic_only_native_directed_inverse_is_certificate=True)


def independent_integral():
    c=MPIntervalContext();c.dps=180;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    owner=current.WholeZSignedLoopFunctions.__new__(current.WholeZSignedLoopFunctions)
    owner.c=c;owner.identity={'manufactured':True};owner.N=160
    owner.owner=lambda ends:SimpleNamespace(flow=f,c=c,P0=[])
    # Independent affine-Z constant source checks the exact kernel/Jacobian.
    rows=dict(kernels={name:f.scalar('-.37') for name in current.RATES},
        Z_derivatives={name:f.scalar('.021') for name in current.RATES})
    owner.query=lambda *args:dict(actual_signed_nonlinear_density_C0_Z=rows)
    comparisons=0
    for chart,left,right in (('O2_axial',(0,1),(1,100)),('O2_slope',(1,4),(1,2)),('O3_power',(0,1),(1,1))):
        got=owner.local_integral(('0','0'),chart,left,right)
        W=p.expm1(p.mpf('.4')) if chart=='O2_axial' else p.mpf('1/4') if chart=='O2_slope' else p.mpf(1)
        lo,hi=current.ep(got['actual_physical_log_width']);assert lo<=W<=hi
        for name,rate in current.RATES.items():
            lam=p.mpf(rate.numerator)/rate.denominator;mass=-p.expm1(-lam*W)/lam if lam else W
            for row,want in zip(got['actual_signed_local_integral_C0_Z'][name],(p.mpf('-.37')*mass,p.mpf('.021')*mass)):
                lo,hi=current.ep(row.finite_interval(max_log=2000));assert lo<=want<=hi;comparisons+=1
        assert not got['actual_full_prefix_C1_defect_functions_installed'] and not got['actual_terminal_controls_installed']
    return dict(passed=True,independent_signed_C0_nonzero_Z_actual_kernel_integrals=comparisons,
        original_axial_exp40_Jacobian_slope_and_power_measures_retained=True,
        local_integrals_do_not_replace_actual_prefix_incoming=True)


def verify(owner,live,counts):
    op=owner.owner(live['exact_Z_cell']);f=op.flow;geometry=live['actual_source_geometry']
    assert live['candidate_N']==owner.N and live['source_identity']==owner.identity
    assert live['exact_common_P0_axial5'] is op.P0
    assert geometry['actual_same_source_Rm_radius'] is op.Rm_factor
    assert geometry['actual_phase_Z_exact_zero'] and geometry['global_phase_not_restarted']
    assert geometry['phase_not_supplied_as_free_angle'] and geometry['adaptive_analytic_interval_digits']>1200
    assert live['original_inverse_functions_evaluated_not_all_u_support_caps']
    assert live['conditional_branches_hulled_after_nonlinear_evaluation_not_added']
    assert live['actual_spatial_Z_is_phase_held_Z_by_radius_Z_zero'] and live['no_y_correction_jet_exported']
    assert live['source_intervals_not_selected_field_values']
    assert all(live[key] is False for key in current.OPEN)
    assert set(live['actual_signed_primitive_C0_Z_phi'])==set(current.OUTPUTS)
    assert not any('_y' in name for name in live['actual_signed_primitive_C0_Z_phi'])
    for alternative in live['actual_signed_inverse_function_alternatives']:
        proof=alternative['original_inverse_and_Z_function_proof']
        assert proof['status']=='enclosed' and proof['original_Z_only_body_and_derivative_algebra_unchanged']
        assert proof['source_roots_and_q_Z_not_conditioned_or_selected']
        if proof['geometry']!='flat':
            inverse=proof['original_C0_inverse'];selected=inverse['selected_inverse']
            assert inverse['inverse_installed_on_this_box'] and not inverse['source_midpoint_used_as_field_value']
            il,ih=current.ep(selected['phase_image']);pl,ph=current.ep(alternative['actual_phase_box'])
            assert il<=ph and ih>=pl
            assert proof['original_direction_is_cosine_Poisson_shape']
            counts['nonflat_directed_inverse_alternatives']+=1
        else:counts['flat_inverse_alternatives']+=1
        for row in alternative['actual_signed_primitive_C0_Z_phi'].values():
            assert row.ctx is f.c and row.scale.bases is f.logs and row.ledger is f.ledger
        density=alternative['actual_nonlinear_density_C0_Z']
        for part in ('kernels','Z_derivatives'):
            current.current.reference.parameters.same_source(f,list(density[part].values()))
        counts['direct_signed_primitive_rows']+=6;counts['nonlinear_signed_density_rows']+=10
    counts['actual_source_function_queries']+=1


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    independent=independent_signed_inverse();integral=independent_integral();radius=scalar_check.radius_identities()
    counts=dict(actual_source_function_queries=0,actual_local_C1_integral_queries=0,
        nonflat_directed_inverse_alternatives=0,flat_inverse_alternatives=0,
        direct_signed_primitive_rows=0,nonlinear_signed_density_rows=0,
        actual_local_integral_C0_Z_rows=0,actual_full_period_integral_phase_covers=0,
        exact_directed_interval_rows=0,typed_rejections=0)
    with mp.workdps(540),patch.object(current.current.WholeZO3RcFiniteN,'contribution',forbidden), \
            patch.object(current.current.original,'OriginalO3RcFiniteN',forbidden), \
            patch.object(current.phase.first,'NativePhaseFirstJets',forbidden), \
            patch.object(current.signed,'NativeSignedUDensityCover',forbidden):
        owner=current.WholeZSignedLoopFunctions()
        assert saved['source_family']==owner.identity and saved['candidate_N']==owner.N==2**3981
        for row in saved['actual_signed_inverse_function_queries']:
            live=owner.query(tuple(row['exact_Z_cell']),row['actual_chart'],tuple(row['exact_left']),tuple(row['exact_right']))
            same_source(encoded(live),row,counts);verify(owner,live,counts)
            assert not live['actual_source_geometry']['full_period']
            for box in live['actual_source_geometry']['actual_phase_boxes']:
                lo,hi=current.ep(box);assert 0<=lo<=hi<=1 and hi-lo<mp.mpf('1e-70')
            print('Current whole-Z signed inverse function audit: '+str(row['exact_Z_cell'])+' '+row['actual_chart'],flush=True)
        for row in saved['actual_signed_local_integral_function_queries']:
            live=owner.local_integral(tuple(row['exact_Z_cell']),row['actual_chart'],tuple(row['exact_left']),tuple(row['exact_right']))
            same_source(encoded(live),row,counts);verify(owner,live['actual_source_function'],counts)
            assert live['live_inverse_function_extension_under_Z_independent_integral']
            assert live['directed_rectangle_function_enclosure_not_selected_integral_value']
            assert live['incoming_correction_not_supplied_or_reset'] and live['physical_Jacobian_applied_once']
            assert not live['actual_full_prefix_C1_defect_functions_installed'] and not live['actual_terminal_controls_installed']
            assert live['actual_source_function']['actual_source_geometry']['full_period']
            op=owner.owner(live['exact_Z_cell']);f=op.flow
            for pair in live['actual_signed_local_integral_C0_Z'].values():
                current.current.reference.parameters.same_source(f,pair)
            if live['actual_chart']=='O3_power':
                assert all(row.zero for pair in live['actual_signed_local_integral_C0_Z'].values() for row in pair)
            counts['actual_local_C1_integral_queries']+=1;counts['actual_local_integral_C0_Z_rows']+=10
            counts['actual_full_period_integral_phase_covers']+=1
        assert counts['actual_source_function_queries']==32 and counts['actual_local_C1_integral_queries']==8
        ends=current.CELLS[0]
        for call in (lambda:owner.radius_query(ends,'O2_slope',(-1,1)),
                lambda:owner.radius_query(ends,'O2_axial',(2,1)),lambda:owner.radius_query(ends,'O2_buffer',(12,1)),
                lambda:owner.radius_query(ends,'O3_power',(3,1)),lambda:owner.radius_query(ends,'O2_slope',(1,1),(0,1)),
                lambda:owner.query(ends,'unknown',(0,1)),lambda:current.exact((1,0)),
                lambda:owner.local_integral(ends,'O3_power',(1,1),(1,1))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid original source function coordinate/integral admitted')
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        genuine_current_whole_Z_signed_inverse_functions_Z_phi_and_nonlinear_densities_checked=True,
        actual_adaptive_global_radius_phase_and_conditional_branch_union_checked=True,
        local_signed_C1_integral_functions_with_true_measure_checked=True,
        independent_signed_inverse_and_Z_diagnostics=independent,independent_actual_integral_measure=integral,
        independent_actual_Rm_radius_identity=radius,replay_counts=counts,
        actual_full_prefix_C1_defect_functions_installed=False,actual_terminal_controls_installed=False,
        actual_global_frequency_admitted=False,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(scalar_check.__file__).name:current.sha(Path(scalar_check.__file__).name),
            Path(scalar_check.original.__file__).name:current.sha(Path(scalar_check.original.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Current whole-Z signed inverse functions and local C1 integral checks passed',flush=True);return receipt


if __name__=='__main__':run()
