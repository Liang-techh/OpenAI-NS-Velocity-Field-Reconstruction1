"""Uniform error bounds for fixed-beta pressure components of the stored datum.

The flatten stage is deliberately excluded. Finite masses are the adapter's
stored zero-center coefficients, interpreted as exact finite numbers.
Runtime evaluation roundoff is not covered by this functional comparison.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def read_interval(ctx, record):
    return ctx.mpf([mp.make_mpf(tuple(record['lower_exact_mpf_tuple'])),
                    mp.make_mpf(tuple(record['upper_exact_mpf_tuple']))])


def fixed_beta_bounds(ctx, mass_error, beta, radius):
    a=ctx.mpf(radius)
    if endpoints(a)[0]<=0 or endpoints(a)[1]>=1:
        raise ValueError('Require 0 < radius < 1')
    if beta not in (0,2):
        raise ValueError('This bound supports only beta 0 or 2')
    lo,hi=endpoints(mass_error)
    size=ctx.mpf(max(abs(lo),abs(hi)))
    factors=[ctx.mpf(1),4*a,4+24*a*a] if beta==2 else [ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
    return [size*f for f in factors]


def run():
    profile=source_profile();e=ScheduleEndpointEnclosures(profile.schedule);iv=e.iv
    receipt=json.loads(Path(__file__).with_name('lei_ren_part1_paper_pointwise_preheat_error_budget.json').read_text())
    adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(e.precision+40):
        finite=adapter.taylor_components(0,center='0')['components']
        check=adapter.taylor_components(2,center=receipt['Z'])['components']
        for stage,row in receipt['stages'].items():
            for n,name in enumerate(('value','first_Z','second_Z')):
                old=mp.make_mpf(tuple(row['finite_positive_integral_derivatives'][name]['exact_mpf_tuple']))
                assert mp.factorial(n)*check[stage]['integral_coefficients'][n]==old, 'Receipt source mismatch: '+stage
        z=iv.mpf(receipt['Z']);g=(1+z*z)**-2
        beta_two={'reference_extension','slope_transition_ref','slope_transition_mu','axial_turnoff','power_buffer','pulse_reserved'}
        names=('value','first_Z','second_Z');rows={};total=[iv.mpf(0) for _ in names]
        for stage,row in receipt['stages'].items():
            if stage=='z_flatten':continue
            beta=2 if stage in beta_two else 0
            truth=read_interval(iv,row['true_positive_integral_intervals']['value'])
            mass=truth/g if beta==2 else truth
            finite_mass=finite[stage]['integral_coefficients'][0]
            error=mass-iv.mpf(finite_mass)
            bounds=fixed_beta_bounds(iv,error,beta,'.8')
            total=[x+y for x,y in zip(total,bounds)]
            rows[stage]=dict(beta=beta,true_mass_interval=mass,finite_mass=finite_mass,
                mass_error_interval=error,uniform_derivative_error_upper_bounds=dict(zip(names,[endpoints(x)[1] for x in bounds])))
        assert len(rows)==13
        report=dict(radius='.8',stage_count=13,stages=rows,
            normalized_uniform_derivative_error_upper_bounds=dict(zip(names,[endpoints(x)[1] for x in total])),
            normalized_uniform_weighted_C2_error_upper=endpoints(total[0]+total[1]+total[2]/2)[1],
            uniform_Z_error_enclosed_for_fixed_beta_stages=True,flatten_included=False,
            full_pressure_error_enclosed=False,finite_datum='stored zero-center coefficient times exact q power',
            adapter_evaluation_roundoff_enclosed=False,original_parameter_errors_enclosed=False,
            core_source_error_enclosed=False,five_defect_interval_closure=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        flatten=finite['z_flatten']
        finite_mass=iv.mpf(0)
        for atom in flatten['atoms']:
            assert 0<=atom['beta']<=2 and atom['atom_at_Z0']>=0
            finite_mass+=iv.mpf(atom['atom_at_Z0'])
        true_mass_upper=e.stage_pressure_upper('z_flatten',axial_radius='.8')['mass_upper']
        # Triangle inequality for positive mixtures; deliberately no relative
        # flatten quadrature accuracy claim is inferred from this coarse bound.
        a=iv.mpf('.8');bound_mass=iv.mpf(true_mass_upper)+finite_mass
        flatten_bounds=[bound_mass,bound_mass*4*a,bound_mass*(4+24*a*a)]
        full=[x+y for x,y in zip(total,flatten_bounds)]
        complete=dict(radius='.8',stage_count=14,
            fixed_beta_receipt=Path(__file__).with_suffix('.json').name,
            flatten_true_mass_upper=true_mass_upper,flatten_finite_mass_interval=finite_mass,
            flatten_uniform_absolute_error_upper_bounds=dict(zip(names,[endpoints(x)[1] for x in flatten_bounds])),
            normalized_uniform_derivative_error_upper_bounds=dict(zip(names,[endpoints(x)[1] for x in full])),
            normalized_uniform_weighted_C2_error_upper=endpoints(full[0]+full[1]+full[2]/2)[1],
            all_pressure_stages_included=True,uniform_Z_error_enclosed=True,
            flatten_method='positive true-plus-finite mass bound; not relative quadrature closure',
            original_parameter_errors_enclosed=False,adapter_evaluation_roundoff_enclosed=False,
            core_source_error_enclosed=False,five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_name('lei_ren_part1_paper_uniform_preheat_error_budget.json').write_text(json.dumps(encode(complete),indent=2)+'\n')
        print('all14 uniform pressure error upper',mp.nstr(endpoints(full[0]+full[1]+full[2]/2)[1],16),flush=True)
        print('13 fixed-beta uniform weighted error upper',mp.nstr(endpoints(total[0]+total[1]+total[2]/2)[1],16),flush=True)


if __name__=='__main__':run()
