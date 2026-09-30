"""Actual Taylor-remainder bounds against retained finite Gauss stage masses."""
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
        for stage in ('slope_transition_ref','slope_transition_mu'):
            coarse=integrator.integrate_stage(stage,panels=32,order=12)
            fine=integrator.integrate_stage(stage,panels=64,order=12)
            assert fine['interval_width']<coarse['interval_width']
            mass=finite['components'][stage]['value_at_Z0'];error=finite_mass_error_interval(fine,mass,e)
            lo,hi=endpoints(error);upper=max(abs(lo),abs(hi))
            a=e.iv.mpf('.8');b=e.iv.mpf(2)
            factors=[e.iv.mpf(1),2*b*a,2*b+4*b*(b+1)*a*a]
            budgets=[endpoints(e.iv.mpf(upper)*f)[1] for f in factors]
            # The exact finite Gauss value is a comparison, not a premise of
            # the analytic integral enclosure.
            rows[stage]=dict(coarse=coarse,fine=fine,finite_Gauss_mass=mass,
                signed_mass_error_interval=error,absolute_mass_error_upper=upper)
            rows[stage]['normalized_axial_error_upper_bounds']=budgets
            rows[stage]['constant_beta']=2
            print(stage,'high-order mass error upper',mp.nstr(upper,16),flush=True)
        report=dict(stages=rows,finite_Gauss_order=192,Taylor_order=12,
            parameter_scope='stored Decimal schedule parameters',
            first_two_stage_mass_errors_enclosed=True,full_pressure_error_enclosed=False,
            core_source_error_enclosed=False,five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':run()
