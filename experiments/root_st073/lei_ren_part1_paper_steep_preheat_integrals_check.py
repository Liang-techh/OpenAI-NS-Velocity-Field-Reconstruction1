"""Actual stored-source steep transitions with high-order integral error bounds."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_variable_preheat_interval_integrals import finite_mass_error_interval
from lei_ren_part1_paper_high_order_preheat_integrals import HighOrderPreheatIntegrals


def run():
    profile=source_profile();e=ScheduleEndpointEnclosures(profile.schedule)
    adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(e.precision+40):
        finite=adapter.taylor_components(0,center=0);integrator=HighOrderPreheatIntegrals(e);rows={}
        for stage in ('steep_transition_in','steep_transition_out'):
            coarse=integrator.integrate_stage(stage,panels=32,order=12)
            fine=integrator.integrate_stage(stage,panels=64,order=12)
            assert fine['interval_width']<coarse['interval_width']
            mass=finite['components'][stage]['value_at_Z0'];error=finite_mass_error_interval(fine,mass,e)
            lo,hi=endpoints(error);upper=max(abs(lo),abs(hi))
            ml,mh=endpoints(fine['normalized_mass_at_Z0_interval'])
            relative=upper/ml
            rows[stage]=dict(coarse=coarse,fine=fine,finite_Gauss_mass=mass,
                signed_mass_error_interval=error,absolute_mass_error_upper=upper,
                relative_mass_error_upper=relative,constant_beta=0,
                normalized_axial_error_upper_bounds=[upper,mp.mpf(0),mp.mpf(0)])
            print(stage,'relative mass error upper',mp.nstr(relative,16),flush=True)
        report=dict(stages=rows,finite_Gauss_order=192,Taylor_order=12,
            parameter_scope='stored Decimal schedule parameters',
            steep_stage_mass_errors_enclosed=True,post_flatten_axial_derivatives_zero=True,
            shared_normalization_preserved=True,full_pressure_error_enclosed=False,
            core_source_error_enclosed=False,five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':run()
