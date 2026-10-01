"""Independent quadrature of constant-power fields and six moments."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_exit_switch_enclosure import constant_power, zero_axial_source_cone
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    c=MPIntervalContext();c.dps=90
    with mp.workdps(130):
        z=IntervalTaylor.variable(c,c.mpf('.3'),1)
        state=dict(phi=1+z/10,U=2+z/5)
        keys=('theta','z','theta_z','p','u_squared','weighted_phi_squared')
        state.update({key:z*0 for key in keys})
        result=constant_power(c,state,c.mpf(3),c.mpf('.08'))
        L=mp.mpf('.08');Z=mp.mpf('.3')
        def reference(key,Z):
            f=1+Z/10;u=2+Z/5
            if key=='phi':return f*mp.exp(-mp.mpf('.4')*L)
            if key=='U':return u
            def integrand(x):
                F=f*mp.exp(-mp.mpf('.4')*x);s=3*mp.exp(x)
                return dict(theta=2*F*s*s,z=u*s,theta_z=2*F*u*s*s,
                    p=F*F*s,u_squared=u*u*s,weighted_phi_squared=F*F*s*s)[key]
            return mp.quad(integrand,[0,L])
        checked=0
        for key in ('phi','U')+keys:
            refs=(reference(key,Z),mp.diff(lambda zz:reference(key,zz),Z))
            for order,v in enumerate(refs):
                lo,hi=endpoints(result[key][order])
                if not lo<=v<=hi:raise AssertionError((key,order,'quadrature outside enclosure'))
                checked+=1
        cutoff=0
        hb=c.mpf('.005')
        for q in ('0','.1','.25','.5','.75','.9','1'):
            qi=c.mpf(q)
            box=1-alpha_box(c,(qi+1)*hb,hb)
            ref=_sigma_mp(mp.mpf(q))
            lo,hi=endpoints(box)
            if not lo<=ref<=hi:raise AssertionError(('cutoff convention',q))
            cutoff+=1
        # Extremely wide nonzero angular shear: algebraic b=0 cancellation
    # must preserve kappa=a and the relaxed margin, not divide repeated boxes.
    class Fixture: pass
    calc=Fixture();calc.ctx=c
    test=zero_axial_source_cone(calc,
        dict(normalized_stress=dict(T_theta_over_F=c.mpf([3,4]),S_z_over_F=c.mpf(0))),
        c.mpf(['1e-1000','.8']))
    if not test['relaxed_cone_certified'] or test['admissible_cone_certified']:
        raise AssertionError('zero-axial-source cone cancellation failed')
    here=Path(__file__).parent
    report=dict(passed=True,fixture_only=True,independent_quadrature_coefficients_contained=checked,
        cutoff_convention_checks=cutoff,zero_axial_shear_cancellation_fixture=True,switch_ODE_independent_numerical_validation=False,
        input_hashes={n:hashlib.sha256((here/n).read_bytes()).hexdigest() for n in
            ('lei_ren_part1_paper_interval_exit_switch_enclosure.py',Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Constant-power quadrature:16 coefficients; cutoff convention:7 points passed',flush=True)
    return report

if __name__=='__main__':run()
