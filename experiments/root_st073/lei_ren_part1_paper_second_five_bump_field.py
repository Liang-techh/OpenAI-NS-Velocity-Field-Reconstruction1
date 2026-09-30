"""Same-source C2 reference moments plus the finite analytic bump inverse."""
import mpmath as mp
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_five_bump_field import _gamma_values,_gamma_derivatives
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


def evaluate_second_bump_field(moment_map,inverse,source,*,Rm,Z,x,delta):
    """Return finite C2 fields, labeled five moments, pressure, Ur and Ur_Z.

Second-Z radial velocity is unavailable without third-Z moment data. No
source at another Z is inferred; caller supplies the same source's C2 data.
"""
    with mp.workdps(moment_map.precision):
        z=mp.mpf(Z);xx=mp.mpf(x);rm=mp.mpf(Rm);dt=mp.mpf(delta)
        if not abs(z)<1 or not 1<=xx<=mp.e or rm<=0 or not 0<=dt<1:
            raise ValueError('Require |Z|<1, 1<=x<=e, Rm>0 and 0<=delta<1')
        if mp.mpf(source['Z'])!=z:raise ValueError('Source axial coordinate differs')
        amplitude=AxialSecondJet(source['Am'],source['Am_Z'],source['Am_ZZ'])
        orders=amplitude.orders
        def constant(v,tangent=0):
            return AxialSecondJet(v,tangent,0,pressure_order=orders[0],width_order=orders[1])
        zd=constant(z,1);R=rm*xx;root=mp.sqrt(2*R);zero=amplitude*0
        scale=mp.sqrt(2)*rm**mp.mpf('1.5')*amplitude
        def physical(rows):
            a1,a2,a3,a4,a5=(rows.get(i,zero) for i in range(1,6))
            return dict(z=rm*a1,theta=scale*a3,theta_z=scale*(a2+4*zd*a3),
                        z_theta=rm*amplitude**2*a4+8*zd*rm*a1,p=amplitude**2*a5)
        base_u=amplitude*xx**mp.mpf('.1')
        parts=dict(reference_power=dict(z=4*zd*R,theta=mp.mpf(5)/8*root*R*base_u,
              theta_z=4*zd*(mp.mpf(5)/8*root*R*base_u),
              z_theta=16*zd*zd*R-mp.mpf(5)/12*R*base_u**2,p=mp.mpf('2.5')*base_u**2))
        for label,rows in source['parts_second_jet'].items():
            if label in parts:raise ValueError('Source part label collision')
            if any(not isinstance(v,AxialSecondJet) for v in rows.values()):
                raise TypeError('C2 field requires complete second-Z source parts')
            parts[label]=physical(rows)
        h=inverse['h'];rows=moment_map.partial(h,amplitude,1,min(xx,mp.mpf(2)))
        parts['five_bump_correction']=physical({i+1:v for i,v in enumerate(rows)})
        moments={name:sum((p[name] for p in parts.values()),zero) for name in parts['reference_power']}
        gammas=_gamma_values(moment_map,xx);primes=_gamma_derivatives(moment_map,xx)
        f=sum((h[i+2]*gammas[i] for i in range(3)),zero)
        fx=sum((h[i+2]*primes[i] for i in range(3)),zero)
        g=h[0]*gammas[0]+h[1]*gammas[2];gx=h[0]*primes[0]+h[1]*primes[2]
        u=amplitude*(xx**mp.mpf('.1')+f);uz=4*zd+g
        uy=amplitude*(mp.mpf('.1')*xx**mp.mpf('.1')+xx*fx);uzy=xx*gx
        P0=source['P0_dual']
        if not isinstance(P0,AxialSecondJet):raise TypeError('Complete common C2 pressure datum required')
        pressure=P0+moments['p'];F=u/root
        m=moments['z'];L=1-dt*z*z
        N=2*z*R*uz.value-(1-dt)*z*m.value-(1-z*z)*m.tangent
        Nz=2*R*(uz.value+z*uz.tangent)-(1-dt)*(m.value+z*m.tangent)+2*z*m.tangent-(1-z*z)*m.second
        ur=N/(L*root);ur_z=Nz/(L*root)+2*dt*z*N/(L*L*root)
        converter=lambda v:PressureWidthJet(v,pressure_order=orders[0],width_order=orders[1])
        stress=evaluate_mp_stress(mp.log(R),z,dt,Utheta=u.value,Uz=uz.value,
            Utheta_y=uy.value,Utheta_Z=u.tangent,Uz_y=uzy.value,Uz_Z=uz.tangent,
            moments={k:v.value for k,v in moments.items()},moments_Z={k:v.tangent for k,v in moments.items()},
            P=pressure.value,P_Z=pressure.tangent,precision=moment_map.precision,
            radius_override=R,scalar_converter=converter,axial_override=converter(z),
            shear_theta=(2*uy.value-u.value)/root,shear_z=root*uzy.value/R)
        return dict(x=xx,Z=z,R=R,second_jet_fields=dict(F=F,Utheta=u,Uz=uz,P=pressure,P0=P0,
                          Utheta_y=uy,Uz_y=uzy),second_jet_moments=moments,
                    second_jet_moment_parts=parts,Ur=ur,Ur_Z=ur_z,Ur_ZZ=None,stress=stress,
                    inverse_increments=inverse['increments'],inverse_terminal_residual=inverse['terminal_residual'],
                    source_part_count=len(source['parts_second_jet']),P0_second_preserved=True,
                    first_Z_radial_velocity_available=True,second_Z_radial_velocity_available=False,
                    uniform_Z_certified=False,functional_closure=False,
                    source_remainder_enclosed=False,infinite_inverse_error_enclosed=False,
                    temporal_recursion=False)
