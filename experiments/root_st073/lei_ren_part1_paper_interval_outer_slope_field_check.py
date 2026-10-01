"""Independent numerical integration and exact endpoint checks for O.2.

The SciPy fixture checks formulas; it is not the directed integral proof.
Production bounds come from monotone rectangles in the field module.
"""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from scipy.integrate import solve_ivp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition,transition_integrals,reference_moments
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    c=MPIntervalContext();c.dps=80;count=0
    def sigma(s):
        if s<=0:return 0.
        if s>=1:return 1.
        a=-1/s**2;b=-1/(1-s)**2;d=b-a
        if d>700:return 0.
        if d < -700:return 1.
        return 1/(1+math.exp(d))
    def ode(s,v):
        j=v[0]
        return [sigma(s),math.exp(1.6*s-.6*j),math.exp(.2*s-1.2*j),math.exp(1.2*s-1.2*j)]
    sol=solve_ivp(ode,(0,1),[0.,0.,0.,0.],method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True)
    if not sol.success:raise AssertionError(sol.message)
    with mp.workdps(90):
        z=IntervalTaylor(c,[c.mpf('.5'),c.mpf(1)])
        A=(1+z*z).reciprocal()*2;P0=z*z+3;Rref=c.mpf(7)
        av=mp.mpf('1.6');az=mp.mpf('-1.28');Z=mp.mpf('.5')
        def contains(box,v):
            nonlocal count
            lo,hi=endpoints(box)
            if not lo<=v<=hi:raise AssertionError(('independent fixture excluded',v,lo,hi))
            count+=1
        for y in ('.25','.5','.75','1'):
            p=evaluate_transition(c,z,c.mpf('.01'),Rref,A,P0,y,256)
            J,Itheta,Ip,Ienergy=map(lambda x:mp.mpf(str(x)),sol.sol(float(y)))
            for bound,value in zip(p['dimensionless_increment_integrals'],(Itheta,Ip,Ienergy)):
                contains(bound,value)
            if y!='1':contains(p['J'],J)
            factor=mp.exp(mp.mpf(y)/10-mp.mpf('.6')*J)
            R=7*mp.exp(mp.mpf(y));u=av*factor;uz=az*factor
            thfactor=mp.sqrt(2)*mp.mpf(7)**mp.mpf('1.5')*(mp.mpf(5)/8+Itheta)
            pfactor=mp.mpf('2.5')+Ip/2
            efactor=7*(mp.mpf(5)/12+Ienergy/2)
            expected=dict(z=(4*Z*R,4*R),theta=(av*thfactor,az*thfactor),
                theta_z=(4*Z*av*thfactor,4*(av+Z*az)*thfactor),
                z_theta=(16*Z*Z*R-av*av*efactor,32*Z*R-2*av*az*efactor),
                p=(av*av*pfactor,2*av*az*pfactor))
            for name,values in expected.items():
                for k,v in enumerate(values):contains(p['physical_moments'][name][k],v)
            # Numerical endpoint field has no rectangle width after exact J(1).
            if y!='1':contains(p['Utheta'][0],u);contains(p['Utheta'][1],uz)
            for k in (0,1):
                expectedP=expected['p'][k]+(mp.mpf('3.25') if k==0 else 1)
                contains(p['P'][k],expectedP)
        inlet=evaluate_transition(c,z,c.mpf('.01'),Rref,A,P0,0)
        exact=reference_moments(c,z,Rref,A)
        for key in exact:
            for k in (0,1):
                if inlet['physical_moments'][key][k]._mpi_!=exact[key][k]._mpi_:
                    # Different grouping can cause a final rounding unit.
                    a,b=endpoints(inlet['physical_moments'][key][k]);d,e=endpoints(exact[key][k])
                    if max(a,d)>min(b,e):raise AssertionError(('inlet discontinuity',key,k))
        endpoint=evaluate_transition(c,z,c.mpf('.01'),Rref,A,P0,1)
        contains(endpoint['J'],mp.mpf('.5'))
        if endpoint['normalized_stress']['S_theta_over_F']._mpi_!=c.mpf(-2)._mpi_:
            raise AssertionError('exact terminal shear lost')
        contains(endpoint['Utheta'][0],av*mp.exp(mp.mpf('-.2')))
        coarse=transition_integrals(c,1,64)[1];fine=transition_integrals(c,1,256)[1]
        for a,b in zip(coarse,fine):
            al,ah=endpoints(a);bl,bh=endpoints(b)
            if not al<=bl<=bh<=ah:raise AssertionError('nested refinement lost')
    here=Path(__file__).parent
    result=dict(passed=True,independent_fixture_containments=count,
        exact_inlet_five_moment_C1_overlap=True,exact_terminal_J_and_shear=True,
        nested_rectangle_refinement=True,
        numerical_oracle='SciPy DOP853 independent coupled primitive ODE, moderate synthetic fixture only',
        numerical_oracle_is_not_a_rigorous_certificate=True,
        input_hashes={n:hashlib.sha256((here/n).read_bytes()).hexdigest() for n in
            (Path(__file__).name,'lei_ren_part1_paper_interval_outer_slope_field.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Outer slope checks passed:',count,'fixture containments; exact endpoints; nested refinement',flush=True)
    return result

if __name__=='__main__':run()
