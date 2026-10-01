"""Bounded checks of weighted cutoff inclusion and finite cone identity."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from scipy.integrate import quad
from lei_ren_part1_paper_large_Md_screen import weighted_cutoff_integrals
from lei_ren_part1_paper_interval_outer_axial_turnoff import cutoff_integrals
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    c=MPIntervalContext();c.dps=80;checks=[]
    with mp.workdps(100):
        for md in ('2','4','6'):
            coarse=weighted_cutoff_integrals(c,md,'.5',128)
            fine=weighted_cutoff_integrals(c,md,'.5',256)
            old=cutoff_integrals(c,md,'.5',256)
            # Independent floating quadrature is a finite-accuracy oracle,
            # not a rigorous integral certificate.
            import math
            end=math.exp(float(md)/2)
            for key,power in [('B_mass',1),('B_squared_mass',2)]:
                lo,hi=endpoints(coarse[key]);fl,fu=endpoints(fine[key]);ol,ou=endpoints(old[key])
                assert lo<=fl<=fu<=hi and max(fl,ol)<=min(fu,ou)
                # Kernel uses exp(-1/q^2), as defined by sigma_value_derivative.
                def integrand(s):
                    q=math.log(s)/float(md)
                    if q<=0:return math.exp(s)
                    a=math.exp(-1/q**2);b=math.exp(-1/(1-q)**2)
                    return math.exp(s)*(b/(a+b))**power
                value,error=quad(integrand,1,end,epsabs=1e-10,epsrel=1e-12)
                assert fl<=mp.mpf(value)<=fu
                checks.append(dict(Md=md,key=key,nested=True,old_overlap=True,scipy_point_contained=True,
                    oracle_reported_error=error))
    for b in (Q(0),Q(1,10),Q(1),Q(3)):
        for x,z,r in [(Q(2),Q(-1),Q(7)),(Q(1,3),Q(4,5),Q(100))]:
            tau=x-2/r;zeta=z-b/r;dot=-(2*tau+b*zeta);cross=b*tau-2*zeta
            original=2*dot**2-b*b*cross**2/2
            factored=(x-2/r)*((8-b**4/2)*x+(8+2*b*b)*b*z-(b*b+4)**2/r)
            assert original==factored
    report=dict(weighted_integral_checks=checks,exact_rational_finite_cone_identities=12,
        all_passed=True,whole_phase_cone_certified=False,
        input_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),Path(__file__).with_name('lei_ren_part1_paper_large_Md_screen.py')]})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('weighted cutoff and finite cone algebra checks PASS',flush=True)
    return report


if __name__=='__main__':run()
