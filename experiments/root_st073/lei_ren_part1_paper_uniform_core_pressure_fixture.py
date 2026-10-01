"""Actual-parameter local perturbations inside the uniform core error bounds."""
import json, math
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets
from lei_ren_part1_paper_core_ra_experiment import _axis_F0_taylor,_axis_u0_taylor
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_component_pressure_core import evaluate_component_core_coefficients
from lei_ren_part1_paper_pressure_core_interval_propagation import core_moment_second_derivatives


def scalar(record):return mp.make_mpf(tuple(record['exact_mpf_tuple']))


def run():
    base=Path(__file__).parent
    bounds=json.loads((base/'lei_ren_part1_paper_uniform_core_pressure_propagation.json').read_text())
    pressure=json.loads((base/'lei_ren_part1_paper_uniform_pressure_high_derivatives.json').read_text())
    precision=473;degree=18;length=21;count=0
    profile=source_profile();datum=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(precision):
        axis=RegularCoreAxisJets(j='1e-14',Lambda='1e36',logC='5e151',delta='1e-200',precision=precision)
        r=mp.mpf(4)/axis.Lambda;physical=mp.exp(28)
        perturb=[physical*scalar(row['normalized_Taylor_coefficient_error_upper'])/2*(-1)**k for k,row in enumerate(pressure['derivatives'][:length])]
        for text in ('-.8','0','.3','.8'):
            z=mp.mpf(text);jets=datum.taylor_components(length-1,center=z)
            nominal=[physical*v for v in jets['pressure_coefficients']]
            kwargs=dict(F0_Z_taylor=_axis_F0_taylor(axis,z,length),U0_Z_taylor=_axis_u0_taylor(z,j=axis.j,length=length),radial_degree=degree,precision=precision)
            a=core_coefficients(z,axis.delta,P0_Z_taylor=nominal,**kwargs)
            b=core_coefficients(z,axis.delta,P0_Z_taylor=[x+y for x,y in zip(nominal,perturb)],**kwargs)
            fa=evaluate_component_core_coefficients(a,r,z,axis.delta)
            fb=evaluate_component_core_coefficients(b,r,z,axis.delta)
            sa=core_moment_second_derivatives(a,r);sb=core_moment_second_derivatives(b,r)
            for name in ('F','Uz','P','Ur'):
                difference=fb[name]-fa[name];upper=scalar(bounds['field_pressure_error_bounds'][name]['absolute_error_upper'])
                assert abs(difference)<=upper;count+=1
            for name in sa:
                for label,key in (('value','moments'),('first_Z','moments_Z'),('second_Z',None)):
                    diff=sb[name]-sa[name] if key is None else fb[key][name]-fa[key][name]
                    upper=scalar(bounds['five_core_moment_pressure_error_bounds'][name][label]['absolute_error_upper'])
                    assert abs(diff)<=upper,(text,name,label);count+=1
    out=dict(all_checks_passed=True,actual_parameter_local_perturbation_checks=count,
        radial_degree=degree,Z_samples=['-.8','0','.3','.8'],point_samples_are_regression_only=True,
        uniform_proof='analytic axis enclosures and paired directed recurrence',radial_series_remainder_enclosed=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Uniform core perturbation fixture passed',count,'actual-parameter checks')


if __name__=='__main__':run()
