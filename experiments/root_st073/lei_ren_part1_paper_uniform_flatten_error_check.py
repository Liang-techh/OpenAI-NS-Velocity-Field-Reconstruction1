"""Stored-schedule directed uniform flatten coefficient receipt."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_uniform_flatten_error import enclose_uniform_flatten_error
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from mpmath.ctx_iv import MPIntervalContext


def run():
    profile=source_profile();adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    result=enclose_uniform_flatten_error(adapter,coefficient_order=96,panels=512,z_radius='.8')
    with mp.workdps(adapter.precision+30):
        mass_lower=endpoints(result['true_mass_interval'])[0]
        iv=MPIntervalContext();iv.dps=adapter.precision+30
        assert mass_lower>0 and result['directed_interval_arithmetic']
        result['uniform_error_relative_to_Z0_mass_lower']={k:endpoints(iv.mpf(v)/iv.mpf(mass_lower))[1] for k,v in result['uniform_error_upper_bounds'].items()}
        result['relative_flat_closure']=False
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n')
        print('Uniform flatten directed receipt; relative value error upper',mp.nstr(result['uniform_error_relative_to_Z0_mass_lower']['value'],16))


if __name__=='__main__':run()
