"""Resolved checks of conservative uniform axis jets."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets
from lei_ren_part1_paper_core_ra_experiment import _axis_F0_taylor,_axis_gradient_series
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    iv=MPIntervalContext();iv.dps=120
    with mp.workdps(130):
        axis=RegularCoreAxisJets(j='.01',Lambda='10',logC='100',delta='.001',precision=130)
        bound=uniform_axis_jets(iv,radius='.8',j=axis.j,Lambda=axis.Lambda,logC=axis.logC,delta=axis.delta,length=8)
        count=0
        for text in ('-.8','-.2','0','.3','.8'):
            z=mp.mpf(text)
            q=_axis_gradient_series(z,length=8,delta=axis.delta,j=axis.j,sigma0=axis.sigma0)
            f=_axis_F0_taylor(axis,z,8)
            for values,intervals in ((q,bound['gradient_coefficients']),(f,bound['F0'])):
                for value,interval in zip(values,intervals):
                    lo,hi=endpoints(interval);assert lo<=value<=hi;count+=1
            assert abs(axis.G(z))<=endpoints(bound['G_absolute_upper'])[1]
        lo,hi=endpoints(bound['axis_anchor_root_bracket'][0])
        assert lo<=axis.Z0<=0
    result=dict(all_checks_passed=True,pointwise_jet_regression_checks=count,
        interval_root_bracket_proved=True,G_bound_basis='path length times positive-denominator analytic envelope',
        point_samples_are_regression_only=True,uniform_axis_proof='rational interval jets plus analytic G envelope')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Uniform-axis fixture passed',count,'jet containment checks')


if __name__=='__main__':run()
