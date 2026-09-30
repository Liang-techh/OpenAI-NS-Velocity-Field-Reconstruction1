"""Termwise field and full quadratic moments of a finite defect response.

Velocity degree N produces quadratic moment degree 2N. No high-degree flux
terms are discarded; this adapter does not claim terminal closure.
"""
from copy import copy
from collections.abc import Mapping
import mpmath as mp
from lei_ren_part1_paper_axial_dual import AxialDual


def _value(v):return v.value if isinstance(v,AxialDual) else v

def _tangent(v):return v.tangent if isinstance(v,AxialDual) else None


def evaluate_response_correction(response,defects,*,Rm,Z,x,delta='1e-200'):
    """Return separate monomial velocity and cumulative-moment increments.

    Rm and x are fixed when differentiating Z. First-Z radial velocity is
    unavailable without second-Z input; U_r is returned as a value only.
    """
    moment_map=response.moment_map
    with mp.workdps(moment_map.precision):
        radius_scale=mp.mpf(Rm);coordinate=mp.mpf(x);z=mp.mpf(Z);dt=mp.mpf(delta)
        if radius_scale<=0 or not 1<=coordinate<=mp.e or not abs(z)<1 or not 0<=dt<1:
            raise ValueError('Require Rm>0, 1<=x<=e, |Z|<1, 0<=delta<1')
        amplitude=response.Am
        if isinstance(amplitude,AxialDual):
            zd=AxialDual(z,1,pressure_order=amplitude.orders[0],width_order=amplitude.orders[1])
        else:zd=z
        radius=radius_scale*coordinate;root=mp.sqrt(2*radius)
        scale=mp.sqrt(2)*radius_scale**mp.mpf('1.5')*amplitude
        endpoint=min(coordinate,mp.mpf(2))
        # Same declared quadrature as the scalar map, restricted to the prefix.
        partial_map=copy(moment_map)
        linear,weights=moment_map._weights_for_interval(mp.mpf(1),endpoint)
        partial_map.linear_matrix=partial_map.matrix=linear
        partial_map.quadratic_weights=weights
        h=response.coefficients
        normalized={}
        for key,vector in h.items():
            normalized[key]=tuple(sum(vector[j]*linear[i,j] for j in range(5)) for i in range(5))
        # Retain all products through 2N, not only the solved formal degrees.
        for ka,a in h.items():
            for kb,b in h.items():
                key=tuple(i+j for i,j in zip(ka,kb))
                q=partial_map.bilinear(a,b,amplitude)
                old=normalized.get(key,(0,)*5)
                normalized[key]=tuple(v+w for v,w in zip(old,q))
        moment_coefficients={}
        for key,rows in normalized.items():
            a1,a2,a3,a4,a5=rows
            moment_coefficients[key]=(scale*a3,radius_scale*a1,
                scale*(a2+4*zd*a3),radius_scale*amplitude**2*a4+8*zd*radius_scale*a1,
                amplitude**2*a5)
        beta=[moment_map.beta(coordinate,c) for c in moment_map.centers]
        velocity_coefficients={}
        for key,vector in h.items():
            f=sum(vector[2+j]*beta[j] for j in range(3))
            g=vector[0]*beta[0]+vector[1]*beta[2]
            ut=amplitude*f
            velocity_coefficients[key]=(ut,g,ut/root)
        # Evaluation is monomial-by-monomial; no aggregate moment subtraction.
        moments=response._evaluate_coefficients(moment_coefficients,defects)
        velocity=response._evaluate_coefficients(velocity_coefficients,defects)
        radial={}
        for key,rows in moments.items():
            dmz=rows[1];dmz_z=_tangent(dmz)
            if dmz_z is None:continue
            g=velocity[key][1] if key in velocity else dmz*0
            radial[key]=(2*z*radius*_value(g)-(1-dt)*z*_value(dmz)-(1-z*z)*dmz_z)/((1-dt*z*z)*root)
        return dict(x=coordinate,R=radius,velocity_terms=velocity,moment_terms=moments,
            radial_velocity_terms=radial,
            moment_order=('theta','z','theta_z','z_theta','p'),
            velocity_order=('Utheta','Uz','F'),
            metadata=dict(velocity_defect_degree=response.degree,quadratic_moment_defect_degree=2*response.degree,
                all_finite_velocity_quadratic_products_retained=True,axis_P0_changed=False,
                radial_first_Z_available=False,source_scope='finite formal response coefficients and supplied defect variables',
                actual_global_field_installed=False,functional_closure=False,cone_certified=False,
                quadrature_enclosed=False,finite_ring_remainder_enclosed=False))
