"""Actual leading core sources for Lei--Ren (13.16)--(13.19).

This materializes the forcing for the first coefficient equations, not the
coefficient solution. Negative-index profiles vanish; Omega_0 has no axial
viscosity term. The radial pressure source is -Omega_0/(2R); full P1_R also
contains 2F0 F1 and cannot be replaced by this source alone.
"""
import json
import math
from pathlib import Path
import numpy as np
from lei_ren_part1_extended_field import ExtendedPartIBackgroundField


def sources(matched,R,Z,*,radial_step,axial_step):
    R,Z=float(R),float(Z); dr,dz=float(radial_step),float(axial_step)
    if not (0<4*dr<R and R+4*dr<matched.R_core and abs(Z)+4*dz<1 and dz>0):
        raise ValueError('Source jets require an interior regular-core point')
    delta=2*matched.core.reference.h
    core=matched.core
    F=lambda r,z:float(core.F(r,z))
    U=lambda r,z:float(core.U(r,z))
    def V_coefficients(z):
        u=np.asarray(core.grid.interpolate(core.u_coefficients,z),float)
        uz=np.asarray(core.grid.interpolate(core.grid.differentiate(core.u_coefficients),z),float)
        powers=np.arange(1,len(u)+1)
        return np.r_[0.,(2*z*u-(1-delta)*z*u/powers-(1-z*z)*uz/powers)/(1-delta*z*z)]
    def V(r,z):
        return float(np.polynomial.polynomial.polyval(r,V_coefficients(z)))
    def derivative(g,r,z,axis):
        h=dr if axis==0 else dz
        call=(lambda k:g(r+k*h,z)) if axis==0 else (lambda k:g(r,z+k*h))
        return (call(-2)-8*call(-1)+8*call(1)-call(2))/(12*h)
    def axial_operator(g,a,r,z):
        return (a*z*g(r,z)+(1-z*z)*derivative(g,r,z,1)
                -2*z*r*derivative(g,r,z,0))/(1-delta*z*z)
    def axial_twice(g,a):
        inner=lambda r,z:axial_operator(g,a,r,z)
        return axial_operator(inner,a-1+delta,R,Z)
    coefficients=V_coefficients(Z)
    v=V(R,Z); vz=derivative(V,R,Z,1)
    vr=float(np.polynomial.polynomial.polyval(R,np.polynomial.polynomial.polyder(coefficients)))
    vrr=float(np.polynomial.polynomial.polyval(R,np.polynomial.polynomial.polyder(coefficients,2)))
    time_v=((1-delta)*Z*vz/2+R*vr)/(1-delta*Z*Z)
    omega=time_v+v*(vr-v/(2*R))+U(R,Z)*axial_operator(V,0,R,Z)-2*R*vrr
    theta=-math.sqrt(2*R)*axial_twice(F,-2-delta)
    axial=-axial_twice(U,-1-delta)
    return {'Omega_0':omega,'radial_pressure_known_source':-omega/(2*R),
            'angular_axial_viscosity_source':theta,'axial_axial_viscosity_source':axial,
            'radial_V_jets':'analytic derivatives of the same finite core polynomial',
            'radial_step':dr,'axial_step':dz,'coefficient_solution_materialized':False}


def run():
    field=ExtendedPartIBackgroundField()
    receipt=json.loads(Path(__file__).with_name('lei_ren_part1_physical_remainder_checks.json').read_text())
    jet_rows=[sources(field.matched,.002,.15,radial_step=.002*f,axial_step=f)
              for f in (2e-4,1e-4)]
    jet=jet_rows[-1]; comparisons=[]
    for row in receipt['rows']:
        if row['R']!=.002 or row['Z']!=.15: continue
        q=row['tau']/(1-row['Z']**2)
        axial_scale=math.sqrt(field.nu)*q**(-1.5+field.h)
        radial_scale=math.sqrt(field.nu)*q**(-1.5)/math.sqrt(2*row['R'])
        measured=row['measurements'][-1]
        predicted=[jet['angular_axial_viscosity_source']*axial_scale,
                   jet['axial_axial_viscosity_source']*axial_scale]
        observed=measured['axial_viscosity_residual'][1:]
        predicted_radial=jet['Omega_0']*radial_scale
        comparisons.append({'tau':row['tau'],'predicted_angular_axial_viscosity':predicted,
            'measured_angular_axial_viscosity':observed,
            'viscosity_max_relative_difference':max(abs(a-b)/max(abs(b),1e-30) for a,b in zip(predicted,observed)),
            'predicted_radial_without_axial_viscosity':predicted_radial,
            'measured_radial_without_axial_viscosity':measured['E_B_without_axial_viscosity'][0],
            'radial_relative_difference':abs(predicted_radial/measured['E_B_without_axial_viscosity'][0]-1)})
    report={'source':'https://arxiv.org/html/2609.35406v1','source_equations':'13.8,13.16-13.19',
            'candidate':field.metadata(),'R':.002,'Z':.15,'source_jet_refinement':jet_rows,
            'independent_physical_comparisons':comparisons,
            'first_coefficient_solved':False,'flatness_established':False,
            'next_equation':'Solve coupled F1,Uz1,P1 with zero axis data and P1_R=2F0F1-Omega0/(2R); then extend and restore positive-order moments.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'jets':jet_rows,'comparisons':comparisons}),flush=True)
    return report


if __name__=='__main__': run()
