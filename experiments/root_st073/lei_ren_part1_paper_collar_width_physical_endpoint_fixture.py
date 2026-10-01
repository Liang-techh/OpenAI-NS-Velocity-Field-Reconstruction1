"""Independent mixed width/Z derivatives of the radial endpoint flux."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet
from lei_ren_part1_paper_collar_width_physical_endpoint import physical_endpoint


def run():
    with mp.workdps(110):
        ctx=MPIntervalContext();ctx.dps=85;z0=mp.mpf('.3');delta=mp.mpf('.01')
        def jet(f):return IntervalAxialSecondJet(ctx,f(z0),mp.diff(f,z0),mp.diff(f,z0,2),pressure_order=0,width_order=2)
        c=lambda v:(lambda z:mp.mpf(v))
        funcs=dict(F=lambda z:2+z/5,Uz=lambda z:mp.mpf('.4')+z/3,
                   mz=lambda z:mp.mpf('.2')+z*z/7,
                   g1=lambda z:-mp.mpf('.2')+z/8,g2=lambda z:mp.mpf('.03')+z*z/9,
                   u1=lambda z:mp.mpf('.04')-z/6,u2=lambda z:-mp.mpf('.01')+z*z/11,
                   mz1=lambda z:mp.mpf('.5')+z*z/4,mz2=lambda z:mp.mpf('.07')-z/5)
        zero=jet(c(0));R=jet(c(2));F=jet(funcs['F'])
        inlet=dict(F=F,R=R,Uz=jet(funcs['Uz']),P=jet(c(3)),P0=jet(c(3)),
                   z=jet(lambda z:z),delta=ctx.mpf(delta),
                   moments={k:zero for k in ('theta','z','theta_z','z_theta','p')})
        inlet['moments']['z']=jet(funcs['mz'])
        first={k:zero for k in ('g','u','theta','mz','mixed','axial','swirl','p')}
        second=dict(first)
        first.update(g=jet(funcs['g1']),u=jet(funcs['u1']),mz=jet(funcs['mz1'])/2)
        second.update(g=jet(funcs['g2']),u=jet(funcs['u2']),mz=jet(funcs['mz2'])/2)
        result=physical_endpoint(inlet,first,second)
        def independent(w,z):
            r=2*mp.exp(2*w);U=funcs['Uz'](z)+w*funcs['u1'](z)+w*w*funcs['u2'](z)
            m=lambda zz:funcs['mz'](zz)+w*funcs['mz1'](zz)+w*w*funcs['mz2'](zz)
            flux=2*z*r*U-(1-delta)*z*m(z)-(1-z*z)*mp.diff(m,z)
            return flux/((1-delta*z*z)*mp.sqrt(2*r))
        rows=[]
        for dz,name in ((0,'Ur'),(1,'Ur_Z')):
            for dw in range(3):
                oracle=mp.diff(lambda w:mp.diff(lambda z:independent(w,z),z0,dz),0,dw)/mp.factorial(dw)
                interval=result[name].component(0,dw).nominal
                lo,hi=(mp.make_mpf(t) for t in interval._mpi_)
                distance=max(lo-oracle,oracle-hi,mp.mpf(0))
                if distance>mp.mpf('1e-80'):raise AssertionError((name,dw,'mixed derivative mismatch'))
                rows.append(dict(field=name,width_power=dw,oracle=mp.nstr(oracle,40),
                                 numerical_oracle_distance=mp.nstr(distance,15)))
        if result['Ur_ZZ'] is not None:raise AssertionError('Unknown third derivative fabricated')
        if result['P0'] is not inlet['P0']:raise AssertionError('Pressure datum replaced')
        report=dict(all_checks_passed=True,mixed_derivative_checks=rows,
                    P0_identity_preserved=True,unknown_Ur_ZZ_preserved=True,
                    actual_source_used=False,full_ODE_error_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print('mixed radial-flux width/Z checks',len(rows),'passed',flush=True)


if __name__=='__main__':run()
