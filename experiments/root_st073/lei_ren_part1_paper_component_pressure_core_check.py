"""Actual-source local component-core receipt; no global installer."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_component_pressure_core import fixture
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    print('building coherent actual source and complete preheat component core',flush=True)
    bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
        logC='5e151',logPstar='14',delta='1e-200',continuous_pressure=True,
        coherent_waiting=True,complete_preheat_components=True)
    core=bundle['component_pressure_core']
    with mp.workdps(core.precision):
        coefficients=core.coefficients('.3')
        z=mp.mpf('.3');delta=core.axis.delta;L=1-delta*z*z
        p0=coefficients['P'][0]
        expected=((1-z*z)*p0[1].component(1)-2*(1+delta)*z*p0[0].component(1))/(2*L)
        actual=coefficients['Uz'][1][0].component(1)
        relative=abs(actual/expected-1)
        assert expected!=0 and relative<mp.mpf('1e-90'),relative
        print('complete pressure propagated through radial degree 18',flush=True)
        samples=[];divergence_errors=[]
        for location in ('1','4'):
            r=mp.mpf(location)/core.axis.Lambda
            values=core.evaluate(r,'.3')
            for k,v in values['divergence_numerator'].atoms.items():
                scale=abs(values['Uz'].component(k))
                if not scale and v:raise AssertionError('Unexpected unscaled divergence component')
                divergence_errors.append(abs(v)/scale if scale else mp.mpf(0))
            samples.append(dict(Lambda_R=location,components={name:
                {str(k):signed_log(v,60) for k,v in poly.atoms.items()}
                for name,poly in values.items() if hasattr(poly,'atoms')},
                moments={name:{str(k):signed_log(v,60) for k,v in poly.atoms.items()}
                    for name,poly in values['moments'].items()}))
        divergence_max=max(divergence_errors,default=mp.mpf(0))
        assert divergence_max<mp.mpf('1e-240'),divergence_max
        prefix=coefficients['P'][0][0].component(0);post=coefficients['P'][0][0].component(1)
        report=dict(Z='.3',radial_degree=18,precision=core.precision,
            datum_prefix=signed_log(prefix,80),datum_post_Rv=signed_log(post,80),
            nominal_sum_loses_tail=(prefix+post)==prefix,
            axial_first_radial_tail_coefficient=signed_log(actual,80),
            independent_U1_replay_relative=mp.nstr(relative,40),
            maximum_retained_pressure_parameter_power=max(k for name in ('F','Uz','P') for row in coefficients[name] for poly in row for k in poly.atoms),
            samples=samples,component_scaled_divergence_max=mp.nstr(divergence_max,40),fixture=fixture(),complete_preheat_axis_input_retained=True,
            pressure_parameter_order_truncated=False,radial_degree_truncated=True,
            post_stage_aggregation_error_enclosed=False,quadrature_error_enclosed=False,
            global_field_installed=False,terminal_pressure_closed=False,
            finite_energy_certified=False,stress_cone_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('datum_post_Rv','nominal_sum_loses_tail',
        'axial_first_radial_tail_coefficient','maximum_retained_pressure_parameter_power')},indent=2),flush=True)
    return report


if __name__=='__main__':run()
