"""Global z-angular-momentum budget, independent of pressure fitting."""
import json
from pathlib import Path
import numpy as np
from .quadrature import unit_rule
from .constrained_tensor_candidate import TensorCandidate
from .constrained_force import RestrictedForce


def angular_moment(field,time,order):
    n,w=unit_rule(order);r=2*n;z=4*n-2
    rr,zz=np.meshgrid(r,z,indexing='ij')
    points=np.stack((rr,np.zeros_like(rr),zz),axis=-1)
    # x*u_y-y*u_x = r*u_theta at theta=0; cylindrical Jacobian r.
    return float(2*np.pi*8*np.sum(rr**2*field(points,time)[...,1]*w[:,None]*w[None,:]))


def run():
    root=Path('artifacts/constrained/tensor_feasible')
    c=TensorCandidate.load(root/'candidate.json')
    f=RestrictedForce(**json.loads((root/'training.json').read_text())['force'])
    rows=[]
    for order in (24,48,96):
        for t in (.25,.5,.75):
            dt=1e-4
            J=lambda t:angular_moment(c.velocity,t,order)
            if t==.25:derivative=(-3*J(t)+4*J(t+dt)-J(t+2*dt))/(2*dt)
            elif t==.75:derivative=(3*J(t)-4*J(t-dt)+J(t-2*dt))/(2*dt)
            else:derivative=(J(t+dt)-J(t-dt))/(2*dt)
            torque=angular_moment(f,t,order)
            mismatch=derivative-torque
            rows.append({'order':order,'time':t,'angular_momentum':J(t),
                'time_derivative':derivative,'force_torque':torque,'budget_mismatch':mismatch,
                'inferred_residual_L2_lower_bound':abs(mismatch)/np.sqrt(512/3)})
    result={'scope':'smooth compact fields on R3; pressure, convection and viscosity have zero integrated torque; assumes exact divergence-free representation',
        'identity':'d/dt integral(x*u_y-y*u_x) = integral(x*f_y-y*f_x) for an exact solution',
        'bound':'Cauchy-Schwarz over [-2,2]^3: L2(residual)>=abs(moment imbalance)/sqrt(512/3)',
        'qualification':'Numerical estimates, not interval-certified lower bounds; dt=1e-4; vary quadrature only',
        'rows':rows}
    (root/'angular_momentum.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([r for r in rows if r['order']==96],indent=2))

if __name__=='__main__':run()
