"""Actual first two pressure transitions: convergent bounds and finite-mass errors."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_variable_preheat_interval_integrals import VariablePreheatIntervalIntegrals,finite_mass_error_interval


def run():
    profile=source_profile();e=ScheduleEndpointEnclosures(profile.schedule)
    adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(e.precision+30):
        finite=adapter.taylor_components(0,center=0)
        integrator=VariablePreheatIntervalIntegrals(e);rows={}
        for stage in ('slope_transition_ref','slope_transition_mu'):
            coarse=integrator.integrate_stage(stage,panels=64)
            fine=integrator.integrate_stage(stage,panels=128)
            assert fine['interval_width']<coarse['interval_width']
            mass=finite['components'][stage]['value_at_Z0']
            lo,hi=endpoints(fine['normalized_mass_at_Z0_interval'])
            assert lo<=mass<=hi, 'Finite mass outside independent radial enclosure'
            error=finite_mass_error_interval(fine,mass,e)
            el,eh=endpoints(error);abs_error=max(abs(el),abs(eh))
            # beta=2 is constant throughout these two stages. Their mass
            # error multiplies exactly q^-2 and hence gives axial budgets.
            a=e.iv.mpf('.8');b=e.iv.mpf(2)
            factors=[e.iv.mpf(1),2*b*a,2*b+4*b*(b+1)*a*a]
            budgets=[endpoints(e.iv.mpf(abs_error)*f)[1] for f in factors]
            rows[stage]=dict(coarse=coarse,fine=fine,finite_Gauss_mass=mass,
                signed_mass_error_interval=error,absolute_mass_error_upper=abs_error,
                normalized_axial_error_upper_bounds=budgets,
                constant_beta=2,finite_Gauss_mass_inside_interval=True)
            print(stage,'mass error upper',mp.nstr(abs_error,12),flush=True)
        negative_rows={}
        for stage in ('axial_turnoff','power_buffer','pulse_reserved','power_buffer_rel','steep_power','waiting'):
            integral=integrator.integrate_negative_stage(stage)
            mass=finite['components'][stage]['value_at_Z0']
            error=finite_mass_error_interval(integral,mass,e)
            el,eh=endpoints(error)
            negative_rows[stage]=dict(integral=integral,finite_mass=mass,
                signed_mass_error_interval=error,absolute_mass_error_upper=max(abs(el),abs(eh)))
        report=dict(axial_radius='.8',stages=rows,negative_slope_stages=negative_rows,finite_Gauss_order=192,
            source_parameter_scope='stored Decimal schedule parameters',
            first_two_stage_mass_errors_enclosed=True,
            full_pressure_error_enclosed=False,core_source_error_enclosed=False,
            five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':run()
