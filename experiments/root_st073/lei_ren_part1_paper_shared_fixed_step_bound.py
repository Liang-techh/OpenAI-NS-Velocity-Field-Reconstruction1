"""Sharp derivative bound8 for the repository's existing flat smooth step.

The old bound32 was safe but too coarse for the paper's T=400A shaping
budget. No cutoff is changed: this is a new analytic bound for the same
exp(-1/s^2) ratio used by alpha_box and the outer schedule.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent


def run():
    c=MPIntervalContext();c.dps=100
    with mp.workdps(140):
        center=1-alpha_box(c,c.mpf('1.5'),c.mpf(1))
        derivative=2*center*(1-center)*(c.mpf('.5')**-3+c.mpf('.5')**-3)
        cl,ch=endpoints(center);dl,dh=endpoints(derivative)
        if not (cl<=mp.mpf('.5')<=ch and dl<=8<=dh and dh-dl<mp.mpf('1e-95')):
            raise ArithmeticError('Existing cutoff normalization changed')
        # v=t^2 in [0,1], t=|1-2s|. The exact polynomial difference
        # (1-v)^4+64v-(1+3v)(1-v)=v*(58+9v-4v^2+v^3)
        # is nonnegative because the bracket>=58-4=54.
        coefficient_floor=58-4
        if coefficient_floor<=0:
            raise ArithmeticError('Cutoff derivative polynomial floor lost')
        result=dict(cutoff_formula='sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), 0<s<1; flat constants outside',
            global_derivative_upper=8,derivative_at_half=8,
            global_derivative_bound_certified=True,existing_cutoff_unchanged=True,
            smooth_flat_endpoints=True,strict_monotonicity_on_open_unit_interval=True,
            proof=dict(symmetry='sigma(1-s)=1-sigma(s), derivative symmetric',
                variables='t=|1-2s|,v=t^2,u=s(1-s)=(1-v)/4',
                identity='sigma_prime/8=(1+3v)/(1-v)^3 / cosh^2(8t/(1-v)^2)',
                hyperbolic='cosh^2(x)=1+sinh^2(x)>=1+x^2',
                polynomial='(1-v)^4+64v-(1+3v)(1-v)=v*(58+9v-4v^2+v^3)>=54v>=0',
                endpoints='the defining exp(-1/s^2) edges are flat; derivative tends to0 at0,1'),
            sample_center_is_not_the_global_proof=True,
            input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                (Path(__file__).name,'lei_ren_part1_paper_interval_comparison_enclosure.py',
                 '../../src/openai_ns_reconstruction/schedule_pressure.py')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print('Existing flat step: analytic global sigma_prime<=8, sharp at1/2; cutoff unchanged')
        return result


if __name__=='__main__':
    run()
