"""New checks of finite-width normalized own-moment transport bounds."""
from pathlib import Path
import sys
import json
import math
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from openai_ns_reconstruction.original_hydro_source_point import finite_width_comparison, width_errors


def run_fixture():
    c=mp.mp.clone();c.dps=100
    iv=MPIntervalContext();iv.dps=160
    controls=[]
    def check(label,condition):
        if not condition:raise AssertionError(label)
        controls.append({'name':label,'passed':True})
    def exact(x):return {'source_enclosure_exact_mpf_tuples':[list(c.mpf(x)._mpf_)]*2}
    domain={'joint_core_domain':{'Cauchy_Z_radius':exact(1),'Cauchy_rho_radius':'.5',
                                'Xh_norm_majorants':{'Psi':exact(4)},'core_value_factor':4}}
    psi=lambda z:c.exp(z)/4+c.sin(z)/8
    slope=lambda z:(1+z*z)/10
    m0=lambda z:c.cos(z)/3
    a0=lambda z:(1+z**4)/5
    h=c.mpf('1e-12');y=2*h;weight=iv.mpf('.2')
    # An independently explicit comparison field and its exact exponential
    # own-history integrals.  This oracle does not call the production ODE.
    moment=lambda z:c.exp(-y)*m0(z)+psi(z)*(1-c.exp(-y))+slope(z)*(y-1+c.exp(-y))
    square_moment=lambda z:(c.exp(-y)*a0(z)+psi(z)**2*(1-c.exp(-y))
        +2*psi(z)*slope(z)*(y-1+c.exp(-y))
        +slope(z)**2*(y*y-2*y+2-2*c.exp(-y)))
    for z in map(c.mpf,('-.3','0','.3')):
        core={name:[iv.mpf(c.diff(fn,z,k)) for k in range(4)]
              for name,fn in (('Psi',psi),('Mpsi',m0),('Apsi',a0))}
        truth,rho,proof=finite_width_comparison(c,iv,core,domain,weight,iv.mpf(h),iv.log(iv.mpf(h)),
                                               source_proved=True,error_cap='1e-3')
        for name,fn in (('Psi',lambda z:psi(z)+y*slope(z)),('Mpsi',moment),('Apsi',square_moment)):
            for k in range(4):
                actual=c.diff(fn,z,k)
                lo,hi=[c.make_mpf(v) for v in truth[name][k]._mpi_]
                check(f'finite_width_{z}_{name}_Z{k}',lo<=actual<=hi)
                # Check the original hb-times-coefficient bound itself,
                # not just the deliberately coarser representable box.
                bound=proof['ordinary_error_proofs'][name][k]['hb_coefficient_bound_exact_mpf_tuples'][1]
                difference=abs(actual-c.diff({'Psi':psi,'Mpsi':m0,'Apsi':a0}[name],z,k))
                check(f'coefficient_bound_{z}_{name}_Z{k}',difference<=h*c.make_mpf(tuple(bound)))
        lo,hi=[c.make_mpf(v) for v in rho._mpi_]
        check(f'physical_rho_{z}',lo<=4*c.exp(2*h)<=hi)
    for name,call in (
        ('missing_source_proof',lambda:width_errors(c,iv,iv.mpf(1),weight,iv.mpf(-10000),orders=2,source_proved=False)),
        ('inadmissible_width_cap',lambda:width_errors(c,iv,iv.mpf(1),weight,iv.mpf(-1),orders=2,source_proved=True)),
        ('negative_error_cap',lambda:width_errors(c,iv,iv.mpf(1),weight,iv.mpf(-10000),orders=2,source_proved=True,error_cap='-1e-3'))):
        try:call()
        except ValueError:check(name,True)
        else:check(name,False)
    return {'passed':True,'control_count':len(controls),'controls':controls,
            'independent_oracle':'explicit affine-in-y comparison field and exact exponential M/A history integrals; mp.diff',
            'original_inverse_or_NS':False,'accepted_ancestor_replayed':False}


if __name__=='__main__':print(json.dumps(run_fixture(),indent=2))
