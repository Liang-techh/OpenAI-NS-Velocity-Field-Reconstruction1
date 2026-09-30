"""Complete common H=1 pressure norm upper bounds for stored schedule data."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    s=source_profile().schedule;e=ScheduleEndpointEnclosures(s)
    with mp.workdps(e.precision+40):
        stages=[name for name in s._stage_bounds if name not in ('reference','heat_connection','exact_heat')]
        rows={name:e.stage_pressure_upper(name,axial_radius='.8') for name in stages}
        terminal=e.terminal_preheat_bounds();iv=e.iv;a=iv.mpf('.8');beta=iv.mpf(2)
        reference=[iv.mpf('2.5')*f for f in (iv.mpf(1),2*beta*a,2*beta+4*beta*(beta+1)*a*a,
                   12*beta*(beta+1)*a+8*beta*(beta+1)*(beta+2)*a**3)]
        total=[reference[i]+sum((iv.mpf(row['normalized_derivative_upper_bounds'][i]) for row in rows.values()),iv.mpf(0))+
               iv.mpf(terminal['normalized_derivative_upper_bounds'][i]) for i in range(4)]
        C2=total[0]+total[1]+total[2]/2
        # Preserve terminal atoms rather than recover them by subtracting totals.
        assert endpoints(terminal['exterior_mass_interval'])[0]>0
        report=dict(precision=e.precision,axial_radius='.8',pre_collar_stages=rows,terminal=terminal,
            reference_derivative_bounds=reference,normalized_derivative_upper_bounds=[endpoints(v)[1] for v in total],
            normalized_C2_upper_bound=endpoints(C2)[1],Pstar_squared=iv.exp(2*e.scalar(s.logPstar)),
            full_stored_parameter_preheat_norm_bounded=True,
            pressure_approximation_error_enclosed=False,original_parameter_errors_enclosed=False,
            core_source_C2_norm_enclosed=False,five_defect_interval_closure=False,
            actual_heat_velocity_field_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('complete stored-parameter preheat norm bounded; normalized C2 upper',mp.nstr(endpoints(C2)[1],16),flush=True)


if __name__=='__main__':run()
