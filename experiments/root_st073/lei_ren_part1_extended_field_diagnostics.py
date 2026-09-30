"""Physical divergence, imposed scale laws and full energy of extended field."""
import json
import math
from pathlib import Path

import numpy as np

from lei_ren_part1_extended_field import ExtendedPartIBackgroundField
from lei_ren_part1_field_diagnostics import gradient, curl, rule
from lei_ren_part1_heat_moments import heat_tail_moments
from openai_ns_reconstruction.coordinates import similarity_coordinates_from_tau


def global_energy(field,tau,*,nR=10,nZ=10):
    near=min(math.sqrt(field.nu)*tau**(.5-field.h),field.z_inner/2)
    zedges=[-field.z_outer,-field.z_inner,-near,0.,near,field.z_inner,field.z_outer]
    ra,rb=field.joined.collar.collar_inner_radius,field.joined.collar.R_b
    redges={0.,.005,.01,.02,ra,rb}
    redges.update(float(r) for r in np.geomspace(.02,ra,18))
    for a,b in field.matched.bump_supports: redges.update((a,b))
    redges=sorted(redges)
    energy=tail_energy=0.
    for a,b in zip(zedges[:-1],zedges[1:]):
        zs,zw=rule(a,b,nZ)
        for z,wz in zip(zs,zw):
            s=similarity_coordinates_from_tau(0.,z/math.sqrt(field.nu),tau,field.h)
            inner=0.
            for left,right in zip(redges[:-1],redges[1:]):
                rs,rw=rule(left,right,nR)
                for R,wr in zip(rs,rw):
                    v=field.cylindrical_at_chart(float(R),s.eta,s.q,float(z),d=s.d)
                    inner+=field.nu*s.q*wr*float(v@v)
            B,_=field.z_cutoff(float(z))
            tail=heat_tail_moments(rb,s.eta,field.joined.collar.heat)
            outer=B*B*field.nu**2*s.q**(-2*field.h)*(-2*tail['quadratic_tail'])
            energy+=math.pi*wz*(inner+outer)
            tail_energy+=math.pi*wz*outer
    return {'total':energy,'radial_heat_tail':tail_energy,
            'inner_connection_and_collar':energy-tail_energy}


def run(field=None):
    field=ExtendedPartIBackgroundField() if field is None else field
    taus=np.array([2.**(-k) for k in range(3,9)])
    rows=[]
    for tau in taus:
        r,z=field.from_similarity(.002,.15,float(tau))
        v=field.cylindrical_from_tau(r,z,float(tau))
        g=gradient(field,np.array([r,0.,z]),float(tau),math.sqrt(field.nu*tau)*1e-4)
        rows.append({'tau':float(tau),'radial_length':r,'axial_length':abs(z),
                     'axial_to_radial_aspect':abs(z)/r,'swirl':abs(v[1]),
                     'axial_velocity':abs(v[2]),'velocity_amplitude':float(np.linalg.norm(v)),
                     'vorticity_amplitude':float(np.linalg.norm(curl(g))),
                     'axial_vorticity':abs(float(curl(g)[2])),
                     'winding_density':abs(v[1]/(2*math.pi*r*v[2])),
                     'global_energy':global_energy(field,float(tau))})
        print(json.dumps({'tau':float(tau),'energy':rows[-1]['global_energy']['total']}),flush=True)
    keys=('radial_length','axial_length','axial_to_radial_aspect','swirl',
          'axial_velocity','velocity_amplitude','vorticity_amplitude',
          'axial_vorticity','winding_density')
    fits={key:float(np.polyfit(np.log(taus),np.log([r[key] for r in rows]),1)[0]) for key in keys}
    tau=1/32
    # Include axis, regular core, cutoff, anchor, four bumps, collar and heat.
    charts=[(0.,.2),(.002,.15),(.012,-.25),(.04,.31)]
    charts.extend(((a+b)/2,.23) for a,b in field.matched.bump_supports)
    charts.append(((field.joined.collar.collar_inner_radius+field.joined.collar.R_b)/2,-.35))
    charts.append((field.joined.collar.R_b*1.1,.37))
    points=[]
    for R,Z in charts:
        r,z=field.from_similarity(R,Z,tau)
        points.append(np.array([r*math.cos(.39),r*math.sin(.39),z]))
    # Fixed axial cutoff contributes to both streamfunction and velocity.
    for z in (.31,-.39):
        s=similarity_coordinates_from_tau(0.,z/math.sqrt(field.nu),tau,field.h)
        points.append(np.array([math.sqrt(2*field.nu*s.q*.012),0.,z]))
    divergence=[]
    for fraction in (1e-3,5e-4,2.5e-4):
        step=math.sqrt(field.nu*tau)*fraction
        gs=[gradient(field,p,tau,step) for p in points]
        errors=[abs(float(np.trace(g))) for g in gs]
        relative=[e/max(1.,float(np.linalg.norm(g))) for e,g in zip(errors,gs)]
        divergence.append({'step':step,'max_abs':max(errors),
                           'max_relative_to_gradient':max(relative),'per_point_abs':errors})
    refined=global_energy(field,float(taus[-1]),nR=16,nZ=16)
    axis=[]
    for eps in (1e-5,5e-6):
        z=field.from_similarity(0.,.2,tau)[1]
        a=field.velocity_from_tau(eps,0.,z,tau)
        b=field.velocity_from_tau(-eps,0.,z,tau)
        center=field.velocity_from_tau(0.,0.,z,tau)
        axis.append({'epsilon':eps,'transverse_odd_error':float(np.max(abs(a[:2]+b[:2]))),
                     'axial_even_error':abs(float(a[2]-b[2])),
                     'center_limit_error':float(np.linalg.norm((a+b)/2-center))})
    report={'metadata':field.metadata(),'rows':rows,'fitted_tau_exponents':fits,
            'divergence':divergence,'divergence_accepted':bool(divergence[-1]['max_abs']<1e-5
                    and divergence[-1]['max_relative_to_gradient']<1e-7),
            'divergence_points':len(points),'axis_checks':axis,
            'outside_axial_support_max':float(np.max(abs(field.velocity_from_tau(.1,.2,field.z_outer,tau)))),
            'finest_time_refined_global_energy':refined,
            'energy_relative_refinement':abs(refined['total']/rows[-1]['global_energy']['total']-1),
            'uniform_energy_integrability_exponent':4*field.h/(1-2*field.h),
            'geometry_scaling_is_imposed_by_representation':True,
            'whole_stress_cone_validated':False,'scale_recursion_established':False,
            'pde_validated':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('fitted_tau_exponents','divergence','divergence_accepted','energy_relative_refinement')}),flush=True)
    return report


if __name__=='__main__': run()
