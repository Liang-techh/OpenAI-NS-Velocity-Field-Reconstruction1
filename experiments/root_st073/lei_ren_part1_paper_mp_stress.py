"""Arbitrary-exponent Section 3.16--3.18 stress from actual profile moments.

All inputs must come from one candidate. No moment targets are substituted.
This evaluator removes binary64 overflow, not the unresolved stress cone.
"""
import mpmath as mp


def evaluate_mp_stress(logR, Z, delta, *, Utheta, Uz, Utheta_y,
                       Utheta_Z, Uz_y, Uz_Z, moments, moments_Z,
                       P, P_Z, precision=160):
    """y=log R derivatives; return MP inertial, shear and total stress."""
    with mp.workdps(precision):
        R=mp.exp(mp.mpf(str(logR))); z=mp.mpf(str(Z)); dt=mp.mpf(str(delta))
        if not abs(z)<1 or not 0<=dt<1:
            raise ValueError('Require |Z|<1 and 0<=delta<1')
        d=1-z*z; L=1-dt*z*z; root=mp.sqrt(2*R)
        a,b,ay,az,by,bz,p,pz=map(mp.mpf,
            (Utheta,Uz,Utheta_y,Utheta_Z,Uz_y,Uz_Z,P,P_Z))
        keys=('theta','z','theta_z','z_theta','p')
        m={k:mp.mpf(moments[k]) for k in keys}
        mz={k:mp.mpf(moments_Z[k]) for k in keys}
        transport=-R+(1-dt)*z*m['z']+d*mz['z']
        Itheta=a*transport/(L*root)+((1-dt/2)*m['theta']
            -(1-dt)*z*mz['theta']/2-d*mz['theta_z']
            +(2*dt-1)*z*m['theta_z'])/(2*L*R)
        Iz=(transport*b+(1-dt)*(m['z']-z*mz['z'])/2
            +2*dt*z*m['z_theta']-d*mz['z_theta']
            +R*(2*(1+dt)*z*p-d*pz))/(L*root)
        Stheta=(2*ay-a)/root; Sz=root*by/R
        Ur=(2*z*R*b-(1-dt)*z*m['z']-d*mz['z'])/(L*root)
        Ur_R=((1+dt)*z*b+2*z*by-d*bz)/(L*root)-Ur/(2*R)
        Nt=-mp.sqrt(R/2)/L*((1+dt)*a/2+(1-dt)*z*az/2+ay)
        Nt-=mp.sqrt(R/2)/L*(-2*(1+dt)*z*b*a
            +d*(bz*a+b*az)-2*z*(by*a+b*ay))
        Nt-=R*Ur_R*a+Ur*ay+Ur*a
        Nz=-mp.sqrt(R/2)/L*((1+dt)*b/2+(1-dt)*z*bz/2+by)
        Nz-=mp.sqrt(R/2)/L*(-2*(1+dt)*z*(b*b+p)
            +d*(2*b*bz+pz)-2*z*(2*b*by+a*a/2))
        Nz-=R*Ur_R*b+Ur*by+Ur*b/2
        return dict(I_theta=Itheta,I_z=Iz,S_theta=Stheta,S_z=Sz,
                    T_theta=Itheta+Stheta,T_z=Iz+Sz,U_r=Ur,
                    N_theta=Nt,N_z=Nz)
