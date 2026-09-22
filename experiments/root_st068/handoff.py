"""Five radial moments and source-leading stress for a LOCAL core.
The result is an interface to a future annulus, not a matched exterior.
"""
import numpy as np
from core_series import add,shift,deriv,product,integrate_x,evaluate


def moment_polynomials(core):
    H=2*shift(core.F,x=1)
    return dict(M=integrate_x(core.U),I=integrate_x(H),J=integrate_x(product(core.U,H)),
                S=integrate_x(add(product(core.U,core.U),-shift(product(core.F,core.F),x=1))),
                Cp=integrate_x(product(core.F,core.F)))


def traces(core,X,eta):
    X,eta=np.broadcast_arrays(np.asarray(X),np.asarray(eta))
    if np.any(X<=0):raise ValueError('Interface radius must be positive')
    p=core.profiles(X,eta);F,U,Pi,v=[p[k][0,0] for k in ['F','U','P','v']]
    h=core.h;A=core.A;D=core.D;d=1-eta**2;L=1-2*h*eta**2
    polys=moment_polynomials(core)
    m={k:evaluate(c,X,eta) for k,c in polys.items()};me={k:evaluate(deriv(c,eta=1),X,eta) for k,c in polys.items()}
    W=1-2*D*eta*m['M']/X-d*me['M']/X;H=2*X*F
    Qs=-W+((1-h)*m['I']-D*eta*me['I']-d*me['J']+2*(h-D)*eta*m['J'])/(X*H)
    Ns=-W*U+(D*(m['M']-eta*me['M'])+4*h*eta*m['S']-d*me['S'])/X+4*A*eta*Pi-d*p['P'][0,1]
    shear_a=-2*X*p['F'][1,0]/F
    stress_theta=X*F*Qs/L+2*X*p['F'][1,0]
    stress_z=np.sqrt(X/2)*(Ns/L+2*p['U'][1,0])
    return dict(X=X,eta=eta,F=F,U=U,Pi=Pi,V0=X*v,F_X=p['F'][1,0],U_X=p['U'][1,0],
                Pi_eta=p['P'][0,1],moments=m,moment_eta_derivatives=me,
                leading_stress_theta=stress_theta,leading_stress_z=stress_z,
                shear_a=shear_a,axial_shear_b=2*X*p['U'][1,0]/(np.sqrt(2*X)*F),
                required_annulus_axial_velocity_integral=-m['M'])
