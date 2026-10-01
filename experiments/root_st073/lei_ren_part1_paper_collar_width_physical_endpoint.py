"""Convert analytic normalized collar width atoms to physical endpoint data.

"Physical" here means unnormalized profile moments and fields, before the
similarity-to-Cartesian space/time map. The formal width variable has unit
atom; coefficients are divided by h_b^w.
No tiny width is materialized. This is a width-order-two endpoint polynomial,
not a completed spatial collar or a certificate for the omitted orders.
"""
from lei_ren_part1_paper_interval_pressure_width_jet import IntervalPressureWidthJet
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet


def physical_moment_coefficients(inlet,normalized):
    r=inlet['R'];Fa=inlet['F'];r2=r*r
    return dict(theta=Fa*r2*normalized['theta'],z=r*normalized['mz'],
                theta_z=Fa*r2*normalized['mixed'],
                z_theta=r*normalized['axial']-Fa*Fa*r2*normalized['swirl'],
                p=Fa*Fa*r*normalized['p'])


def physical_endpoint(inlet,first,second,*,s=2):
    """Return same-pressure endpoint fields and all five moments through W^2."""
    template=inlet['F'];ctx=template.ctx;orders=dict(pressure_order=template.pressure_order,
                                                  width_order=template.width_order)
    if template.width_order!=2:
        raise ValueError('This endpoint adapter requires width order two')
    W=IntervalAxialSecondJet(ctx,IntervalPressureWidthJet(ctx,{(0,1):1},**orders))
    first_m=physical_moment_coefficients(inlet,first)
    second_m=physical_moment_coefficients(inlet,second)
    moments={name:inlet['moments'][name]+W*first_m[name]+W*W*second_m[name]
             for name in first_m}
    R=inlet['R']*(1+W*s+W*W*(s*s)/2)
    F=inlet['F']*(1+W*first['g']+W*W*(second['g']+first['g']*first['g']/2))
    U=inlet['Uz']+W*first['u']+W*W*second['u']
    P=inlet['P']+W*first_m['p']+W*W*second_m['p']
    # Ur needs the first moment derivative; Ur_Z additionally needs its second.
    # The unknown third derivative is deliberately not promoted to zero.
    z=inlet['z'].value;delta=IntervalPressureWidthJet(ctx,inlet['delta'],**orders)
    r=R.value;root=(2*r).sqrt();d=1-z*z;L=1-delta*z*z
    mz=moments['z'];flux=2*z*r*U.value-(1-delta)*z*mz.value-d*mz.tangent
    flux_Z=2*r*U.value+2*z*r*U.tangent-(1-delta)*mz.value+(1+delta)*z*mz.tangent-d*mz.second
    Ur=flux/(L*root)
    Ur_Z=flux_Z/(L*root)+2*delta*z*flux/(L*L*root)
    return dict(R=R,F=F,Uz=U,P=P,moments=moments,
                first_width_physical_moments=first_m,second_width_physical_moments=second_m,
                Ur=Ur,Ur_Z=Ur_Z,Ur_ZZ=None,P0=inlet['P0'],
                width_atom_convention='coefficient divided by h_b**width_power',
                original_pressure_datum_preserved=True,width_order=2,
                endpoint_only=True,omitted_width_orders_enclosed=False,
                full_spatial_divergence_certified=False,source_errors_enclosed=False,
                temporal_recursion=False)
