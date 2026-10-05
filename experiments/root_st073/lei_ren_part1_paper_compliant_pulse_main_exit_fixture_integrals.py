"""Independent moderate-parameter original gp/beta integrals, computed once.

This numerical fixture is not production data or a cone certificate. A hash
binds its generating source and precision so unchanged quadrature is reused.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp

PATH=Path(__file__)
PARAMETERS=dict(mu='.2',ell='.15',logG='log(.2*.07)',precision=100)


def run():
    source=hashlib.sha256(PATH.read_bytes()).hexdigest()
    output=PATH.with_suffix('.json')
    if output.exists():
        old=json.loads(output.read_bytes())
        if old.get('generator_sha256')==source and old.get('parameters')==PARAMETERS and old.get('mpmath_version')==mp.__version__:
            return old
    with mp.workdps(100):
        mu=mp.mpf(PARAMETERS['mu']);ell=mp.mpf(PARAMETERS['ell']);logG=mp.log(mu*mp.mpf('.07'))
        sigma=lambda x:mp.mpf(0) if x<=0 else mp.mpf(1) if x>=1 else 1/(1+mp.exp(1/x**2-1/(1-x)**2))
        def gp(x):
            if x<=0 or x>=11:return mp.mpf(0)
            primitive=x-mp.mpf('.01') if x>=mp.mpf('.02') else mp.quad(
                lambda v:sigma(50*v),[0,x/2,x],method='gauss-legendre')
            return primitive*sigma(11-x)
        normal=mp.quad(lambda v:mp.exp(-1/(1-v*v)),[-1,0,1])
        beta=lambda v:mp.exp(-1/(1-(v/ell)**2))/(ell*normal) if abs(v)<ell else mp.mpf(0)
        W=[mp.quad(lambda v,lam=mp.mpf('.5')-i*mu:mp.exp(lam*v)*beta(v),[-ell,0,ell]) for i in (1,2)]
        gram=mp.quad(lambda v:mp.exp(-2*mu*v)*beta(v)**2,[-ell,0,ell])
        A=[mp.exp((mp.mpf('.5')-mu)*center)*W[0] for center in (-3,-1)]
        row2=[mp.exp((mp.mpf('.5')-2*mu)*center)*W[1] for center in (-3,-1)]
        D=[(a-b)/mu for a,b in zip(row2,A)]
        weights=[mp.exp(-26-2*mu*center)*gram for center in (-3,-1)]
        entrance_linear=[]
        for i in (1,2):
            k=(mp.mpf('.5')-i*mu)/mu;a=mp.mpf('.02')
            entrance_linear.append((mp.exp(k*a)*mp.mpf('.01')
                -mp.quad(lambda v:mp.exp(k*v)*sigma(50*v),[0,a/2,a]))/k)
        entrance_energy=mp.quad(lambda a:mp.exp(-2*a)*gp(a)**2,
            [0,mp.mpf('.01'),mp.mpf('.02')],method='gauss-legendre')
        energy_primitive=lambda a:-mp.exp(-2*a)*((a-mp.mpf('.01'))**2/2+(a-mp.mpf('.01'))/2+mp.mpf('.25'))
        exit_energy=mp.quad(lambda a:mp.exp(-2*a)*gp(a)**2,[10,mp.mpf('10.5'),11])
        Kpulse=entrance_energy+energy_primitive(10)-energy_primitive(mp.mpf('.02'))+exit_energy
        def linear_integral(x,i):
            k=(mp.mpf('.5')-i*mu)/mu
            prim=lambda a:mp.exp(k*a)*((a-mp.mpf('.01'))/k-1/k**2)
            return (entrance_linear[i-1]+prim(10)-prim(mp.mpf('.02'))
                +mp.quad(lambda a:mp.exp(k*a)*gp(a),[10,(10+x)/2,x]))/mu
        full_rows=[mp.exp(-13*(mp.mpf('.5')-i*mu)/mu)*linear_integral(mp.mpf(11),i)/mp.exp(logG) for i in (1,2)]
        values=dict(normal=normal,W=W,gram=gram,A=A,D=D,weights=weights,entrance_linear=entrance_linear,
                    entrance_energy=entrance_energy,exit_energy=exit_energy,Kpulse=Kpulse,full_rows=full_rows)
        def encode(value):
            if isinstance(value,list):return [encode(v) for v in value]
            return mp.nstr(value,110)
        result=dict(parameters=PARAMETERS,mpmath_version=mp.__version__,generator_sha256=source,
                    values={key:encode(value) for key,value in values.items()},
                    fixture_only=True,production_data=False)
        output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('Independent original gp/beta fixture integrals generated and source-bound',flush=True)
        return result


if __name__=='__main__':run()
