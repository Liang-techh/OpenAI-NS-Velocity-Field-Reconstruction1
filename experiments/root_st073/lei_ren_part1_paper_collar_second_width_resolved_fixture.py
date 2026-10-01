"""Independent resolved collar ODE check of analytic width derivatives.

Uses a smooth switch and prescribed analytic core drivers, not actual source
data. RK convergence is a formula check, not a certified RK remainder.
"""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet
from lei_ren_part1_paper_collar_first_width_interval_endpoint import first_width_coefficients
from lei_ren_part1_paper_collar_second_width_moments import second_width_coefficients


def run():
    with mp.workdps(80):
        R,F,u,D,Iz,DR,IzR,a=map(mp.mpf,('2','3','.4','.7','.2','.03','-.02','.1'))
        sigma_integral=mp.quad(lambda x:x*_sigma_mp(x),[0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1])
        K=mp.mpf('.5')+sigma_integral;iv=MPIntervalContext();iv.dps=75
        def j(x):return IntervalAxialSecondJet(iv,x,pressure_order=0,width_order=0)
        inlet=dict(R=j(R),F=j(F),Uz=j(u),D=j(D),I_z=j(Iz),D_R=j(DR),I_z_R=j(IzR),F_R=j(a*F/R))
        names=('g','u','theta','mz','mixed','axial','swirl','p')
        first=first_width_coefficients(inlet);second=second_width_coefficients(inlet,iv.mpf(K))
        def midpoint(jet):
            pair=jet.value.component(0,0);v=pair.nominal
            return (mp.make_mpf(v._mpi_[0])+mp.make_mpf(v._mpi_[1]))/2
        target1=[midpoint(first[n]) for n in names];target2=[midpoint(second[n]) for n in names]
        epsilon=mp.mpf('1e-18');rows=[]
        for steps in (128,256):
            h=mp.mpf(2)/steps
            # Switch nodes shared by both signed-width resolved integrations.
            switches=[mp.mpf(1) if k==0 else (mp.mpf(0) if k>=steps else _sigma_mp(1-mp.mpf(k)/steps))
                      for k in range(2*steps+1)]
            def integrate(W):
                state=[mp.mpf(0),u]+[mp.mpf(0)]*6
                def rhs(s,y,chi0):
                    g,v=y[:2];radius=R*mp.exp(s*W)
                    chi=chi0+W*(1-chi0);barD=D+DR*(radius-R);barIz=Iz+IzR*(radius-R)
                    return [W*(-chi*barD/2),
                        W*(-chi*mp.sqrt(radius/2)*mp.exp(g-a*s*W)*barIz),
                        W*2*mp.exp(2*s*W+g),W*mp.exp(s*W)*v,
                        W*2*mp.exp(2*s*W+g)*v,W*mp.exp(s*W)*v*v,
                        W*mp.exp(2*s*W+2*g),W*mp.exp(s*W+2*g)]
                for n in range(steps):
                    s=n*h;k1=rhs(s,state,switches[2*n])
                    k2=rhs(s+h/2,[v+h*q/2 for v,q in zip(state,k1)],switches[2*n+1])
                    k3=rhs(s+h/2,[v+h*q/2 for v,q in zip(state,k2)],switches[2*n+1])
                    k4=rhs(s+h,[v+h*q for v,q in zip(state,k3)],switches[2*n+2])
                    state=[v+h*(p+2*q+2*r+t)/6 for v,p,q,r,t in zip(state,k1,k2,k3,k4)]
                return state
            plus=integrate(epsilon);minus=integrate(-epsilon);zero=[mp.mpf(0),u]+[mp.mpf(0)]*6
            c1=[(p-m)/(2*epsilon) for p,m in zip(plus,minus)]
            c2=[(p+m-2*z)/(2*epsilon**2) for p,m,z in zip(plus,minus,zero)]
            rows.append(dict(steps=steps,first_width_errors=[mp.nstr(abs(x-y),30) for x,y in zip(c1,target1)],
                             second_width_errors=[mp.nstr(abs(x-y),30) for x,y in zip(c2,target2)]))
        previous=max(mp.mpf(x) for x in rows[0]['second_width_errors'])
        latest=max(mp.mpf(x) for x in rows[1]['second_width_errors'])
        if not latest<mp.mpf('1e-7') or not latest<previous/8:
            raise AssertionError('Resolved width coefficient convergence failed')
        result=dict(all_checks_passed=True,states=list(names),resolved_signed_width=mp.nstr(epsilon),
            rows=rows,second_width_error_contraction=mp.nstr(previous/latest,20),
            source_data_used=False,certified_RK_remainder=False,full_ODE_error_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
        print('second-width independent resolved error',mp.nstr(latest,12),
              'contraction',mp.nstr(previous/latest,12),flush=True)


if __name__=='__main__':run()
