"""Pressure-independent circulation diagnostic for a fixed velocity and force."""
import json
from pathlib import Path
import numpy as np
from .constrained_validation import residual
from .quadrature import unit_rule


def loop_circulation(velocity,force,rectangle,time=.75,order=64,step=.0025):
    r0,r1,z0,z1=rectangle;n,w=unit_rule(order);r=r0+(r1-r0)*n;z=z0+(z1-z0)*n
    if not (0<r0<r1 and z0<z1):raise ValueError('ordered positive-radius rectangle required')
    edges=[(np.column_stack((r,np.zeros(order),np.full(order,z0))),0,r1-r0),
           (np.column_stack((np.full(order,r1),np.zeros(order),z)),2,z1-z0),
           (np.column_stack((r,np.zeros(order),np.full(order,z1))),0,-(r1-r0)),
           (np.column_stack((np.full(order,r0),np.zeros(order),z)),2,-(z1-z0))]
    total=0.
    for points,axis,length in edges:
        v=residual(velocity,lambda p,t:np.zeros(np.shape(p)[:-1]),force,points,time,step=step)['momentum']
        total+=length*np.sum(w*v[:,axis])
    return float(total)


def run():
    from .constrained_local_pressure import LocalPressureCandidate
    root=Path('artifacts/constrained/local_pressure');c=LocalPressureCandidate.load(root/'candidate.json')
    rectangle=[.6,1.,.2,.6];perimeter=1.6;rows=[]
    for order,h in [(96,.005),(96,.0025),(96,.00125),(32,.00125),(64,.00125)]:
        value=loop_circulation(c.velocity,c.force,rectangle,order=order,step=h)
        rows.append({'quadrature_order':order,'derivative_step':h,'circulation':value,'estimated_sup_residual_lower_bound':abs(value)/perimeter})
    report={'scope':'numerical pressure-independent estimate; not interval-certified','time':.75,
            'rectangle_r0_r1_z0_z1':rectangle,'selection':'largest estimate among25 predefined .4-wide development rectangles',
            'identity':'integral grad(p) dot dl=0 on closed rectangle; sup residual >=abs(circulation)/perimeter',
            'refinement':'first vary derivative step at order96, then quadrature at fixed step.00125','rows':rows}
    (root/'circulation.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
