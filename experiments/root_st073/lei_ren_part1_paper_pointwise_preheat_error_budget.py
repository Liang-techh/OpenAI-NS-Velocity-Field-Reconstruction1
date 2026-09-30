"""All-stage pressure approximation error budget for one common stored source."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_high_order_preheat_integrals import HighOrderPreheatIntegrals
from lei_ren_part1_paper_variable_preheat_interval_integrals import VariablePreheatIntervalIntegrals
from lei_ren_part1_paper_flatten_preheat_integrals import integrate_flatten


def run():
    profile=source_profile();e=ScheduleEndpointEnclosures(profile.schedule);iv=e.iv
    adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(e.precision+40):
        Z='.3';z=iv.mpf(Z);q=1+z*z
        factors=[q**-2,-4*z*q**-3,-4*q**-3+24*z*z*q**-4]
        data=adapter.taylor_components(2,center=Z)
        high=HighOrderPreheatIntegrals(e);negative=VariablePreheatIntervalIntegrals(e)
        names=('value','first_Z','second_Z');rows={}
        def add(stage,true,method):
            finite=[mp.factorial(n)*data['components'][stage]['integral_coefficients'][n] for n in range(3)]
            errors=[truth-iv.mpf(v) for truth,v in zip(true,finite)]
            upper=[max(abs(lo),abs(hi)) for lo,hi in map(endpoints,errors)]
            rows[stage]=dict(true_positive_integral_intervals=dict(zip(names,true)),
                finite_positive_integral_derivatives=dict(zip(names,finite)),
                signed_pressure_error_intervals=dict(zip(names,[-v for v in errors])),
                absolute_error_upper_bounds=dict(zip(names,upper)),method=method)
        add('reference_extension',[iv.mpf('2.5')*v for v in factors],'exact reference integral')
        for stage in ('slope_transition_ref','slope_transition_mu','steep_transition_in','steep_transition_out'):
            mass=high.integrate_stage(stage,panels=64,order=12)['normalized_mass_at_Z0_interval']
            true=[mass*f for f in factors] if stage.startswith('slope_') else [mass,iv.mpf(0),iv.mpf(0)]
            add(stage,true,'order12 interval Taylor radial integral')
        for stage in ('axial_turnoff','power_buffer','pulse_reserved','power_buffer_rel','steep_power','waiting'):
            mass=negative.integrate_negative_stage(stage)['normalized_mass_at_Z0_interval']
            true=[mass*f for f in factors] if stage in ('axial_turnoff','power_buffer','pulse_reserved') else [mass,iv.mpf(0),iv.mpf(0)]
            add(stage,true,'negative-slope exponential lower/upper integral')
        flatten=integrate_flatten(e,Z=Z,panels=256,order=12)
        add('z_flatten',[flatten['integrals'][name] for name in names],'order12 variable-beta integral, pointwise Z')
        terminal=e.terminal_preheat_bounds()
        add('heat_collar',[terminal['heat_collar_mass_interval'],iv.mpf(0),iv.mpf(0)],'H=1 collar K0 interval')
        add('exterior_power_tail',[terminal['exterior_mass_interval'],iv.mpf(0),iv.mpf(0)],'exact H=1 infinite power integral')
        assert set(rows)==set(data['components']), 'A pressure stage is missing'
        total=[sum((iv.mpf(row['absolute_error_upper_bounds'][name]) for row in rows.values()),iv.mpf(0)) for name in names]
        C2=total[0]+total[1]+total[2]/2
        report=dict(Z=Z,precision=e.precision,stages=rows,stage_count=len(rows),
            normalized_pointwise_derivative_error_upper_bounds=dict(zip(names,[endpoints(v)[1] for v in total])),
            normalized_pointwise_weighted_C2_error_upper=endpoints(C2)[1],
            physical_pressure_scale=iv.exp(2*e.scalar(profile.schedule.logPstar)),
            all_pressure_stage_errors_included=True,source_profile='same shared continuous preheat schedule',
            pointwise_Z_only=True,uniform_Z_error_enclosed=False,
            original_parameter_errors_enclosed=False,core_source_error_enclosed=False,
            five_defect_interval_closure=False,actual_heat_velocity_field_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('all-stage pointwise pressure C2 error upper',mp.nstr(endpoints(C2)[1],16),'stages',len(rows),flush=True)


if __name__=='__main__':run()
