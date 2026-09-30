"""Actual flatten mass and two axial derivatives against common Gauss data."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_flatten_preheat_integrals import integrate_flatten


def run():
    profile=source_profile();e=ScheduleEndpointEnclosures(profile.schedule)
    adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(e.precision+40):
        z=mp.mpf('.3');receipt=adapter.taylor_components(2,center=z)
        row=receipt['components']['z_flatten'];finite={name:mp.factorial(n)*row['integral_coefficients'][n]
            for n,name in enumerate(('value','first_Z','second_Z'))}
        coarse=integrate_flatten(e,Z='.3',panels=128,order=12)
        fine=integrate_flatten(e,Z='.3',panels=256,order=12)
        errors={}
        for name in finite:
            assert fine['interval_widths'][name]<coarse['interval_widths'][name]
            interval=fine['integrals'][name]
            error=interval-e.iv.mpf(finite[name]);lo,hi=endpoints(error)
            upper=max(abs(lo),abs(hi));relative=upper/abs(finite[name])
            errors[name]=dict(finite_Gauss_value=finite[name],signed_error_interval=error,
                absolute_error_upper=upper,relative_error_upper=relative)
            print('flatten',name,'relative error upper',mp.nstr(relative,16),flush=True)
        for result in (coarse,fine):
            result['stored_left']=str(result['stored_left']);result['stored_right']=str(result['stored_right'])
        report=dict(Z='.3',coarse=coarse,fine=fine,errors=errors,finite_Gauss_order=192,
            pointwise_value_and_axial_errors_enclosed=True,uniform_Z_errors_enclosed=False,
            original_parameter_errors_enclosed=False,full_pressure_error_enclosed=False,
            core_source_error_enclosed=False,five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':run()
