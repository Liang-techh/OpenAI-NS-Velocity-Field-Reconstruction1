"""Independent moderate-size algebra checks for the physical norm ledger.

Fixtures verify the weighted reciprocal bound and the Cstar powers in
the exact five-moment stress formula. They are not production profiles.
"""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_shared_physical_norm_family import LogBounds
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent


def run():
    with mp.workdps(75):
        c = MPIntervalContext(); c.dps = 65; logs = LogBounds(c)
        sum_checks = 0
        for values in ([mp.mpf('.001'), 1, 100], [2, 3], [mp.mpf('1e-100'), mp.mpf('1e100')]):
            major = logs.add(*(logs.number(v) for v in values))
            if mp.log(sum(values)) > endpoints(major)[1]:
                raise AssertionError('Log sum envelope lost a positive term')
            sum_checks += 1
        reciprocal_checks = 0
        q = lambda s, z: 2+mp.mpf('.3')*s+mp.mpf('.2')*z*z
        # On s in [0,2], Z in [-1,1]: q>=2 and weighted norm<=3.7.
        for m in (3, 4):
            upper = mp.exp(endpoints(logs.reciprocal(logs.number('3.7'), c.mpf(2), m))[1])
            for s in (0, mp.mpf('.7'), 2):
                for z in (-1, mp.mpf('.2'), 1):
                    norm = sum(abs(mp.diff(lambda a, b: 1/q(a,b), (s,z), (i,k)))
                               / (math.factorial(i)*math.factorial(k))
                               for i in range(m+1) for k in range(m+1-i))
                    if norm > upper:
                        raise AssertionError('Reciprocal mixed jet exceeded majorant')
                    reciprocal_checks += 1

        # Polynomial normalized core. r is deliberately ordinary so the
        # original physical formula can be evaluated without cancellation.
        r = mp.mpf('.04'); dt = mp.mpf('.001')
        def moments(R, z, C):
            a = mp.mpf('.8')+mp.mpf('.01')*z*z
            b = mp.mpf('-.1'); v0 = 4*z+mp.mpf('.01'); v1 = mp.mpf('.02')*z
            f = (a+b)/C; v = v0+v1
            theta = 2*r*r*(a/2+b/3)/C
            mz = r*(v0+v1/2)
            tz = 2*r*r*(a*v0/2+(a*v1+b*v0)/3+b*v1/4)/C
            ztheta = r*(v0*v0+v0*v1+v1*v1/3)-r*r*(a*a/2+2*a*b/3+b*b/4)/(C*C)
            pressure = r*(a*a+a*b+b*b/3)/(C*C)
            return dict(f=f, v=v, theta=theta+f*(R*R-r*r),
                z=mz+v*(R-r), theta_z=tz+f*v*(R*R-r*r),
                ztheta=ztheta+v*v*(R-r)-f*f*(R*R-r*r)/2,
                p=pressure+f*f*(R-r))

        def stress(R, z, C):
            m = moments(R,z,C); L = 1-dt*z*z; d = 1-z*z
            dz = lambda key: mp.diff(lambda t: moments(R,t,C)[key],z)
            A = (1-dt)*z*m['z']+d*dz('z')
            angular = m['f']*(-R+A)/L
            angular += ((1-dt/2)*m['theta']-(1-dt)*z*dz('theta')/2
                        -d*dz('theta_z')+(2*dt-1)*z*m['theta_z'])/(2*L*R)
            P0 = lambda t: -30/(1+t*t)**2
            P = P0(z)+m['p']; Pz = mp.diff(P0,z)+dz('p')
            axial = ((-R+A)*m['v']+(1-dt)*(m['z']-z*dz('z'))/2
                     +2*dt*z*m['ztheta']-d*dz('ztheta')
                     +R*(2*(1+dt)*z*P-d*Pz))/(L*mp.sqrt(2*R))
            return angular/m['f'], axial/m['f']

        scaling_checks = 0
        for R in (mp.mpf('.04'), mp.mpf('.3'), mp.mpf(110)):
            for z in (mp.mpf('-.8'), 0, mp.mpf('.6')):
                D2,E2 = stress(R,z,2); D5,E5 = stress(R,z,5)
                # E(C)=C*E0+C^-1*E_swirl for fixed normalized data.
                E0 = (5*E5-2*E2)/21
                Es = 2*E2-4*E0
                D11,E11 = stress(R,z,11)
                scale = max(1,abs(D2),abs(E11))
                if max(abs(D5-D2),abs(D11-D2),abs(E11-11*E0-Es/11)) > mp.mpf('1e-60')*scale:
                    raise AssertionError('Physical Cstar powers do not cancel as claimed')
                scaling_checks += 3
        result = dict(log_sum_checks=sum_checks,reciprocal_mixed_jet_checks=reciprocal_checks,
            physical_five_moment_Cstar_power_checks=scaling_checks,
            moderate_fixture_is_algebra_only=True,production_family_not_sampled=True,
            all_checks_passed=True,full_NS_validation=False,
            input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                (Path(__file__).name,'lei_ren_part1_paper_shared_physical_norm_family.py')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print('PASS: 3 logarithmic sums, 18 mixed reciprocal jets, 27 exact physical Cstar-power checks')
        return result


if __name__ == '__main__':
    run()
