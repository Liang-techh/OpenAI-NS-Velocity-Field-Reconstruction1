"""Accepted coherent source: sharper directed flatten coefficient enclosure."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_affine_flatten_coefficients import enclose_affine_flatten_error
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    profile,alignment=accepted_profile();adapter=ContinuousPreheatPressure(profile,quadrature_order=192)
    result=enclose_affine_flatten_error(adapter,coefficient_order=96,panels=512,z_radius='.8')
    with mp.workdps(300):
        iv=MPIntervalContext();iv.dps=300
        lower=endpoints(result['true_mass_interval'])[0];assert lower>0
        result['uniform_error_relative_to_Z0_mass_lower']={k:endpoints(iv.mpf(v)/iv.mpf(lower))[1] for k,v in result['uniform_error_upper_bounds'].items()}
        result['accepted_schedule_sha256']=alignment['accepted_schedule']['sha256']
        result['relative_flat_closure']=False
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n')
        for name,value in result['uniform_error_relative_to_Z0_mass_lower'].items():print(name,'uniform relative-to-mass error upper',mp.nstr(value,16),flush=True)


if __name__=='__main__':run()
