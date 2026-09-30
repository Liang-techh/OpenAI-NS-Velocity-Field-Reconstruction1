"""Physical divergence, multiscale geometry and full radial-tail energy.

Scale fits are measurements of this fixed self-similar representation, not
evidence of recursive PDE closure. No residual-defined force is introduced.
"""
import json
import math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from lei_ren_part1_field import PartIBackgroundField
from lei_ren_part1_heat_moments import heat_tail_moments
from openai_ns_reconstruction.coordinates import similarity_coordinates_from_tau


def rule(a, b, n):
    x, w = leggauss(n)
    return a+(b-a)*(x+1)/2, w*(b-a)/2


def gradient(field, point, tau, step):
    columns = []
    for axis in range(3):
        off = np.zeros(3); off[axis] = step
        fn = lambda p: field.velocity_from_tau(*p, tau)
        columns.append((fn(point-2*off)-8*fn(point-off)
                        +8*fn(point+off)-fn(point+2*off))/(12*step))
    return np.column_stack(columns)


def curl(g):
    return np.array([g[2,1]-g[1,2], g[0,2]-g[2,0], g[1,0]-g[0,1]])


def sector_energy(field, tau, nR=12, nZ=24):
    rs, rw = rule(0., .05, nR)
    zs, zw = rule(-.5, .5, nZ)
    energy = 0.
    for Z, wz in zip(zs, zw):
        q = tau/(1-Z*Z)
        z = math.sqrt(field.nu)*q**(.5-field.h)*Z
        dz = math.sqrt(field.nu)*q**(.5-field.h)*(1-2*field.h*Z*Z)/(1-Z*Z)
        for R, wr in zip(rs, rw):
            v = field.cylindrical_at_chart(float(R), float(Z), q, z)
            energy += math.pi*field.nu*q*dz*wr*wz*float(v@v)
    return energy


def global_energy(field, tau, nR=12, nZ=14):
    # Split physical z at the contracting core scale and both cutoff edges.
    near = min(math.sqrt(field.nu)*tau**(.5-field.h), field.z_inner/2)
    edges = [-field.z_outer, -field.z_inner, -near, 0., near,
             field.z_inner, field.z_outer]
    radial_edges = [0., .05, .055, .095, .115, .18, .2]
    energy = 0.
    tail_energy = 0.
    for za, zb in zip(edges[:-1], edges[1:]):
        zs, zw = rule(za, zb, nZ)
        for z, wz in zip(zs, zw):
            s = similarity_coordinates_from_tau(0., z/math.sqrt(field.nu), tau, field.h)
            inner = 0.
            for ra, rb in zip(radial_edges[:-1], radial_edges[1:]):
                rs, rw = rule(ra, rb, nR)
                for R, wr in zip(rs, rw):
                    v = field.cylindrical_at_chart(float(R), s.eta, s.q, float(z), d=s.d)
                    inner += field.nu*s.q*wr*float(v@v)
            B, _ = field.z_cutoff(float(z))
            tail = heat_tail_moments(field.joined.R_join, s.eta, field.joined.heat)
            outer = B*B*field.nu**2*s.q**(-2*field.h)*(-2*tail['quadratic_tail'])
            energy += math.pi*wz*(inner+outer)
            tail_energy += math.pi*wz*outer
    return dict(total=energy, radial_heat_tail=tail_energy, matched_inner_and_annulus=energy-tail_energy)


def run():
    field = PartIBackgroundField()
    taus = np.array([2.**(-k) for k in range(3,9)])
    rows = []
    for tau in taus:
        r, z = field.from_similarity(.02, .15, float(tau))
        v = field.cylindrical_from_tau(r, z, float(tau))
        g = gradient(field, np.array([r,0.,z]), float(tau), math.sqrt(field.nu*tau)*1e-4)
        winding = abs(v[1]/(2*math.pi*r*v[2]))
        rows.append(dict(tau=float(tau), radial_length=r, axial_length=abs(z),
                         axial_to_radial_aspect=abs(z)/r, swirl=abs(v[1]), axial_velocity=abs(v[2]),
                         velocity_amplitude=float(np.linalg.norm(v)), vorticity_amplitude=float(np.linalg.norm(curl(g))),
                         axial_vorticity=abs(float(curl(g)[2])), winding_density=winding,
                         winding_over_scaled_axial_length=winding*abs(z),
                         moving_core_sector_energy=sector_energy(field,float(tau)),
                         global_energy=global_energy(field,float(tau))))
    keys = ('radial_length','axial_length','axial_to_radial_aspect','swirl','axial_velocity',
            'velocity_amplitude','axial_vorticity','vorticity_amplitude','winding_density',
            'winding_over_scaled_axial_length','moving_core_sector_energy')
    fits = {key: float(np.polyfit(np.log(taus),np.log([r[key] for r in rows]),1)[0]) for key in keys}
    tau = 1/32
    points = []
    for R,Z in ((.015,.07),(.075,.15),(.13,.31),(.19,-.22),(.25,.37),(0.,.2)):
        r,z=field.from_similarity(R,Z,tau)
        points.append(np.array([r*math.cos(.39),r*math.sin(.39),z]))
    for z in (.31,-.39):
        s=similarity_coordinates_from_tau(0.,z/math.sqrt(field.nu),tau,field.h)
        points.append(np.array([math.sqrt(2*field.nu*s.q*.1),0.,z]))
    divergence = []
    for fraction in (1e-3,5e-4,2.5e-4):
        step=math.sqrt(field.nu*tau)*fraction
        gs=[gradient(field,p,tau,step) for p in points]
        defects=[abs(float(np.trace(g))) for g in gs]
        normalized=[a/max(1.,float(np.linalg.norm(g))) for a,g in zip(defects,gs)]
        divergence.append(dict(step=step,max_abs=max(defects),max_relative_to_gradient=max(normalized),
                               per_point_abs=defects))
    refined_energy=global_energy(field,float(taus[-1]),nR=20,nZ=22)
    coarse_energy=rows[-1]['global_energy']['total']
    axis_checks=[]
    for eps in (1e-5,5e-6):
        z=field.from_similarity(0.,.2,tau)[1]
        a=field.velocity_from_tau(eps,0.,z,tau)
        b=field.velocity_from_tau(-eps,0.,z,tau)
        center=field.velocity_from_tau(0.,0.,z,tau)
        axis_checks.append(dict(epsilon=eps,transverse_odd_error=float(np.max(abs(a[:2]+b[:2]))),
                                axial_even_error=abs(float(a[2]-b[2])),
                                center_limit_error=float(np.linalg.norm((a+b)/2-center))))
    outside=field.velocity_from_tau(.1,.2,field.z_outer,tau)
    report=dict(metadata=field.metadata(), rows=rows, fitted_tau_exponents=fits,
                expected_exponents=dict(radial_length=.5, axial_length=.5-field.h,
                                        axial_to_radial_aspect=-field.h, swirl=-.5-field.h,
                                        axial_velocity=-.5-field.h, axial_vorticity=-1-field.h,
                                        winding_density=-.5, winding_over_scaled_axial_length=-field.h),
                divergence=divergence,
                divergence_accepted=bool(divergence[-1]['max_abs']<1e-5 and
                                         divergence[-1]['max_relative_to_gradient']<1e-7),
                divergence_tolerances=dict(absolute=1e-5,relative_to_gradient=1e-7),
                axis_smoothness_checks=axis_checks,
                axial_support_outside_max=float(np.max(abs(outside))),
                finest_time_refined_global_energy=refined_energy,
                finest_time_energy_relative_refinement=abs(refined_energy['total']/coarse_energy-1),
                uniform_energy_integrability_exponent=4*field.h/(1-2*field.h),
                energy_bound_scope='For fixed h<1/6, bounded profiles and zero axial primitive outside Rjoin, q^(-2h)<=|z/sqrt(nu)|^(-4h/(1-2h)) away from z=0. This integrable bound plus H<=1 controls the full radial heat tail and fixed axial cutoff uniformly as tau decreases. No PDE/stress bound follows.',
                geometry_scaling_is_imposed_by_representation=True,
                pde_validated=False,scale_recursion_established=False,
                stress_resolved=False,full_five_moment_matching=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(fits=fits,divergence=divergence,
                         global_energies=[r['global_energy']['total'] for r in rows],
                         energy_refinement=report['finest_time_energy_relative_refinement'])))
    return report


if __name__=='__main__':
    run()
