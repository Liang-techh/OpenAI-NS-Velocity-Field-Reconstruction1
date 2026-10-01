"""Scalar physical-stress comparison and strict/uncertain cone fixtures."""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_interval_exit_stress_enclosure import evaluate,cone
from lei_ren_part1_paper_interval_comparison_enclosure_check import fixture
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,constant
from lei_ren_part1_paper_candidate_shared_inlet import finite_moments
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    with mp.workdps(120):
        calc = fixture()
        c = calc.ctx
        z = IntervalTaylor.variable(c,c.mpf('.3'),3)
        calc.z,calc.center_family = z,'.3'
        calc.F0 = c.mpf('.7')
        calc.S = constant(c,calc.F0*calc.F0,3)
        psi = [1+z/5,constant(c,'.1',3)]
        urows = [2+z/4,constant(c,'-.1',3)]
        s,R = c.mpf(2),c.mpf('.2')
        phi,u = psi[0]+psi[1]*s,urows[0]+urows[1]*s
        m = finite_moments(psi,urows,s)
        state = dict(phi=phi.truncate(1),U=u.truncate(1),
                     **{key:value.truncate(1) for key,value in m.items()})
        gy = (psi[1]*s/phi).truncate(1)
        uy = (urows[1]*s).truncate(1)
        packet = dict(normalized_exit_state=state,physical_R=R,
                      actual_terminal_g_y=gy,actual_terminal_U_y=uy)
        bounded = evaluate(calc,packet)
        def mid(v):
            lo,hi = endpoints(v)
            return (lo+hi)/2
        F = phi*calc.F0
        physical = dict(theta=m['theta']*calc.F0*calc.eps**2,z=m['z']*calc.eps,
            theta_z=m['theta_z']*calc.F0*calc.eps**2,
            z_theta=m['u_squared']*calc.eps-calc.S*m['weighted_phi_squared']*calc.eps**2,
            p=calc.S*m['p']*calc.eps)
        P = calc.p0+physical['p']
        r,zv = mid(R),mid(z[0])
        root = mp.sqrt(2*r)
        fv,uv,g0,uy0 = mid(F[0]),mid(u[0]),mid(gy[0]),mid(uy[0])
        reference = evaluate_mp_stress(mp.log(r),zv,mid(calc.delta),
            Utheta=root*fv,Uz=uv,Utheta_y=root*fv*(mp.mpf('.5')+g0),
            Utheta_Z=root*mid(F[1]),Uz_y=uy0,Uz_Z=mid(u[1]),
            moments={key:mid(value[0]) for key,value in physical.items()},
            moments_Z={key:mid(value[1]) for key,value in physical.items()},
            P=mid(P[0]),P_Z=mid(P[1]),precision=120,
            shear_theta=2*fv*g0,shear_z=root*uy0/r)
        for key,value in bounded['stress'].items():
            lo,hi = endpoints(value)
            if not lo <= reference[key] <= hi:
                raise AssertionError(('physical stress reference outside interval',key))
        relaxed = cone(c,c.mpf(-1),c.mpf(0),c.mpf(2),c.mpf(0))
        strong = cone(c,c.mpf(-3),c.mpf(0),c.mpf(1),c.mpf(0))
        failed = cone(c,c.mpf(-1),c.mpf(0),c.mpf(-1),c.mpf(0))
        uncertain = cone(c,c.mpf(-3),c.mpf(0),c.mpf([-1,1]),c.mpf(0))
        branch = cone(c,c.mpf([-3,-1]),c.mpf(0),c.mpf(1),c.mpf(0))
        if not relaxed['relaxed_cone_certified'] or relaxed['admissible_cone_certified']:
            raise AssertionError('Relaxed cone fixture misclassified')
        if not strong['admissible_cone_certified'] or not failed['family_cone_ruled_out']:
            raise AssertionError('Strict cone fixture misclassified')
        if uncertain['relaxed_cone_certified'] or uncertain['family_cone_ruled_out'] or branch['branch']!='not isolated':
            raise AssertionError('Uncertain interval mistaken for certificate or failure')
        out = dict(passed=True,fixture_only=True,scalar_physical_stress_components_contained=6,
            cone_fixture_cases=5,uncertain_intervals_not_claimed_failed=True,
            input_hashes={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                for name in ('lei_ren_part1_paper_interval_exit_stress_enclosure.py',
                             'lei_ren_part1_paper_interval_exit_stress_enclosure_check.py')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
        print('Directed stress fixtures:6 physical components,5 cone classifications passed',flush=True)
        return out


if __name__ == '__main__':
    run()
