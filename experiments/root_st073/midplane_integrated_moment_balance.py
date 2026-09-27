"""Axisymmetric momentum moments using integrated conservation identities."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from adaptive_bridge_recursive_defect import build_fields
from joined_field import coordinates
from midplane_resolved_feasibility import RADIAL_BREAKS
from midplane_outer_moment_dae import LocalTimeCorrection
from separated_moment_modes import SeparatedMomentModes,RADIAL_WINDOWS_THREE
from radial_continuation import ROOT


def integral_state(field,inner,R,z,tau,order):
    q=float(coordinates(0.,z/np.sqrt(inner.nu),tau,inner.h)['q'])
    ri=np.sqrt(2*inner.nu*q*inner.p.X_max)
    edges=sorted(set(np.clip([0.,ri,R,*[ri*(1+15*f) for f in RADIAL_BREAKS]],0.,R)))
    g,w=leggauss(order);r=np.concatenate([(a+b)/2+(b-a)*g/2 for a,b in zip(edges[:-1],edges[1:])])
    weights=np.concatenate([(b-a)*w/2 for a,b in zip(edges[:-1],edges[1:])])
    points=np.column_stack((r,np.zeros_like(r),np.full_like(r,z)))
    u,p=field.fields(points,tau)
    return np.array([weights@(r*r*u[:,1]),weights@(r*u[:,2]),
        weights@(r*r*u[:,2]*u[:,1]),weights@(r*u[:,2]**2),weights@(r*p)])


def balance(field,inner,R,z,tau,hz,ht,hr,order):
    center=integral_state(field,inner,R,z,tau,order)
    zm2,zm,zp,zp2=[integral_state(field,inner,R,z+j*hz,tau,order) for j in (-2,-1,1,2)]
    tm2,tm,tp,tp2=[integral_state(field,inner,R,z,tau+j*ht,order) for j in (-2,-1,1,2)]
    dz=(zm2-8*zm+8*zp-zp2)/(12*hz)
    dzz=(-zp2+16*zp-30*center+16*zm-zm2)/(12*hz*hz)
    dtau=(tm2-8*tm+8*tp-tp2)/(12*ht)
    pts=np.array([[R+j*hr,0.,z] for j in (-2,-1,0,1,2)])
    u,_=field.fields(pts,tau);ur,ut,uz=u[2]
    dr=(u[0]-8*u[1]+8*u[3]-u[4])/(12*hr)
    theta=-dtau[0]+dz[2]+R*R*ur*ut-field.nu*(R*R*dr[1]-R*ut+dzz[0])
    axial=-dtau[1]+dz[3]+R*ur*uz+dz[4]-field.nu*(R*dr[2]+dzz[1])
    return np.array([-theta/(R*R),-axial/R])


def evaluate(field,inner,k,order,zfactor=.002):
    tau=.5*2.**-k;out=[]
    for eta in (-.2,.2):
        pt=inner.from_similarity([inner.p.X_max*16**2],[eta],tau)[0];R,z=float(pt[0]),float(pt[2])
        zp=inner.from_similarity([inner.p.X_max],[eta+.01],tau)[0,2]
        zm=inner.from_similarity([inner.p.X_max],[eta-.01],tau)[0,2]
        hz=zfactor*abs(zp-zm)/.02
        out.extend(balance(field,inner,R,z,tau,hz,1e-4*tau,5e-4*np.sqrt(inner.nu*tau),order))
    return np.array(out)


def run():
    inner,fields=build_fields();base=fields['two_sided_cone']
    class Polynomial:
        nu=inner.nu
        def fields(self,pts,tau):
            x,y,z=np.asarray(pts).T;a=.3;b=.7;c=.2
            return np.column_stack((a*x-b*y,a*y+b*x,-2*a*z)),c*z*z
    check=balance(Polynomial(),inner,2.,.3,.7,.002,.0001,.0001,24)
    expected=np.array([-.3*.7*2**2/2,-(2*.3**2+.2)*.3*2])
    np.testing.assert_allclose(check,expected,rtol=1e-7,atol=1e-8)
    class TimePolynomial:
        nu=inner.nu
        def fields(self,pts,tau):
            x,y,z=np.asarray(pts).T;a=.3+.1*tau;b=.7-.2*tau
            return np.column_stack((a*x-b*y,a*y+b*x,-2*a*z)),.2*z*z
    time_check=balance(TimePolynomial(),inner,2.,.3,.7,.002,.0001,.0001,24)
    aa=.3+.1*.7;bb=.7-.2*.7
    time_expected=np.array([-(.2+2*aa*bb),-(.1+2*aa*aa+.2)*.3*2])
    np.testing.assert_allclose(time_check,time_expected,rtol=1e-7,atol=1e-8)
    print('Static and time-dependent polynomial identity checks passed',flush=True)
    a=np.array(json.loads((ROOT/'midplane_outer_axial_interscale.json').read_text())['amplitudes'])
    current=SeparatedMomentModes(base,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    report=dict(polynomial_check=dict(computed=check.tolist(),expected=expected.tolist()),time_polynomial_check=dict(computed=time_check.tolist(),expected=time_expected.tolist()),scales=[],
        scope='Axisymmetric integrated conservation identity with fixed integration radius under differentiation; regular-axis boundary terms assumed. Numerical derivatives of integrated quantities remain to be controlled.',accepted=False,scale_recursion_established=False)
    for k in (13,17):
        source=json.loads((ROOT/f'midplane_outer_moment_dae_k{k}.json').read_text())
        f=LocalTimeCorrection(current,source['coefficient_values'],source['coefficient_slopes'],k)
        rows=[]
        for order,zfactor in ((48,.002),(96,.002),(96,.001)):
            m=evaluate(f,inner,k,order,zfactor)
            row=dict(order=order,z_step_factor=zfactor,moments=m.tolist(),absolute_max=float(max(abs(m))))
            rows.append(row);print(json.dumps(dict(k=k,**row)),flush=True)
        report['scales'].append(dict(k=k,replay=rows,previous_cartesian_replay=source['replay']))
        (ROOT/'midplane_integrated_moment_balance.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())

if __name__=='__main__':
    run()
