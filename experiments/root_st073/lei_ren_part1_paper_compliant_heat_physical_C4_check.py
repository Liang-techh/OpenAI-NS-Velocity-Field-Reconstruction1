"""Independent Cartesian momentum and viscosity-pullback exterior checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_heat_physical_C4 import CompliantHeatPhysicalC4, physical_heat_identities, viscosity_source_row
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'heat_physical_C4.json'


def independent_Cartesian_momentum_fixture():
    """Full Gamma and pressure integral, genuine Cartesian derivatives.

    Moderate amplitudes exercise operators; source transfer is proved
    separately. No actual extreme amplitude or global residual is fitted.
    """
    with mp.workdps(45):
        a=mp.mpf('.15'); delta=2*a; A=mp.mpf('1.3')
        residuals=[];divergences=[]
        for nu,tau,x,y in ((mp.mpf('.01'),mp.mpf('.7'),mp.mpf('1.2'),mp.mpf('.9')),
                           (mp.mpf('.7'),mp.mpf('.4'),mp.mpf('.8'),mp.mpf('1.1'))):
            amplitude=A*nu**(1+a)
            def H(xi):
                if not xi:return mp.mpf(1)
                return xi**(-1-a)*mp.hyperu(1+a,2,1/xi)
            def velocity(component,xx,yy,zz,remaining):
                r=mp.sqrt(xx*xx+yy*yy)
                swirl=amplitude*r**(-1-delta)*H(4*nu*remaining/r**2)
                return (-yy*swirl/r if component==0 else xx*swirl/r if component==1 else mp.mpf(0))
            point=(x,y,mp.mpf('.3')); vals=[velocity(i,*point,tau) for i in range(3)]
            r=mp.sqrt(x*x+y*y);xi=4*nu*tau/r**2
            def Hp(v):
                if not v:return -a*(1+a)
                return -a*(1+a)*v**(-2-a)*mp.hyperu(2+a,2,1/v)
            # Differentiate the entire pressure integral under the integral
            # sign once. Nested numerical differentiation of quadrature
            # unnecessarily doubles precision and repeats this same tail.
            B=mp.quad(lambda q:q**delta*H(xi*q)**2,[0,'.25',1])
            Bprime=2*mp.quad(lambda q:q**(delta+1)*H(xi*q)*Hp(xi*q),[0,'.25',1])
            pressure_r=amplitude**2*r**(-3-2*delta)*((1+delta)*B+xi*Bprime)
            divergence=mp.mpf(0)
            for component in range(3):
                gradients=[];laplace=mp.mpf(0)
                for axis in range(3):
                    def fn(q):
                        p=list(point);p[axis]=q
                        return velocity(component,*p,tau)
                    gradients.append(mp.diff(fn,point[axis]))
                    laplace+=mp.diff(fn,point[axis],2)
                divergence+=gradients[component]
                pr=pressure_r*point[component]/r if component<2 else mp.mpf(0)
                dt=-mp.diff(lambda remaining:velocity(component,*point,remaining),tau)
                residual=dt+sum(vals[i]*gradients[i] for i in range(3))+pr-nu*laplace
                if abs(residual)>mp.mpf('1e-34'):raise ArithmeticError('Independent Cartesian exterior momentum failed')
                residuals.append(mp.nstr(abs(residual),12))
            if abs(divergence)>mp.mpf('1e-34'):raise ArithmeticError('Independent Cartesian exterior divergence failed')
            divergences.append(mp.nstr(abs(divergence),12))
        return dict(physical_viscosities=['.01','.7'],Cartesian_momentum_components_checked=6,
                    full_infinite_pressure_integral_differentiated=True,
                    Cartesian_divergences_checked=2, absolute_momentum_residuals=residuals,
                    absolute_divergence_residuals=divergences,moderate_fixture_only=True,passed=True)


def independent_viscosity_row_fixture():
    c=MPIntervalContext(); c.dps=100
    nu=c.mpf('.01'); lognu=c.ln(nu); count=0
    base=dict(exact_zero=False,log_absolute_upper=c.ln(3),terms=[dict(log_absolute_upper=c.ln(3),signed_coefficient=c.mpf(3))])
    for amplitude_power in (1,2):
        for degree in range(5):
            power=c.mpf(amplitude_power-degree)/2
            row=viscosity_source_row(base,lognu,power)
            target=c.ln(3*nu**power)
            lo,hi=endpoints(row['log_absolute_upper']);tl,th=endpoints(target)
            if max(lo,tl)>min(hi,th) or row['physical_viscosity_exponent']!=power:
                raise ArithmeticError('Physical viscosity derivative units lost')
            count+=1
    return dict(spatial_viscosity_unit_rows_checked=count,passed=True)


def run():
    record=json.loads((HERE/NAME).read_bytes()); hashes=dict(record['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Physical exterior source changed: '+name)
    field=CompliantHeatPhysicalC4()
    if field.bridge!=record['physical_source_bridge'] or physical_heat_identities()!=record['physical_heat_identities']:
        raise ValueError('Actual physical source/identity bridge changed')
    c=MPIntervalContext();c.dps=240;read=lambda value:read_interval(c,value)
    zeros=finite=0
    for point in record['samples']+[record['whole_source_exterior']]:
        for gate in ('actual_compliant_field_physical_heat_region_certified','regional_physical_viscosity_leading_NS_identity_certified',
                     'viscosity_source_pullback_verified','original_pressure_datum_and_positive_amplitude_retained'):
            if not point[gate]:raise ValueError('Actual physical exterior gate missing: '+gate)
        if point['actual_physical_map_and_prefactor_transfer_pending']:raise ValueError('Physical source transfer still pending')
        for value in list(point['completed_background_stress_tensor_components'].values())+list(point['physical_momentum_residual_components'].values())+[
                point['physical_axial_viscosity_remainder'],point['physical_radial_remainder'],point['physical_divergence']]:
            if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Actual regional exact zero identity lost')
            zeros+=1
        for component,parts in point['zeroth_physical_source_rows'].items():
            for row in parts.values():
                if not row['exact_zero'] and not mp.isfinite(endpoints(read(row['log_absolute_upper']))[1]):
                    raise ArithmeticError('Nonfinite physical heat source bound')
                finite+=1
        for flag in ('global_admissible_stress_lift_constructed','whole_outer_cone_certified','independently_bounded_global_flat_remainder',
                     'full_background_NS_validation','physical_energy_integral_certified','temporal_recursion'):
            if point[flag]:raise ValueError('Regional heat identity promoted globally: '+flag)
    # Exercise an actual non-unit viscosity source packet, not only a
    # moderate synthetic derivative fixture.
    with mp.workdps(280):
        actual=field.exterior('.5','4',log_tau='-10',theta='.7',viscosity='.01')
        # Compare the requested decimal above the interval context's
        # precision, rather than a default53-bit rounded representative.
        if not endpoints(actual['physical_viscosity'])[0]<=mp.mpf('.01')<=endpoints(actual['physical_viscosity'])[1]:
            raise ValueError('Physical viscosity selection changed')
    fixture=independent_Cartesian_momentum_fixture(); units=independent_viscosity_row_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],implicit_source_sha256=record['implicit_source_sha256'],
                actual_compliant_field_physical_heat_region_certified=True,actual_physical_map_and_prefactor_transfer_pending=False,
                regional_physical_viscosity_leading_NS_identity_certified=True,
                actual_source_physical_zero_bounds_checked=zeros,actual_finite_zeroth_source_rows_checked=finite,
                actual_nonunit_viscosity_packet_checked=True,physical_source_bridge=field.bridge,
                physical_heat_identities=field.proof,independent_Cartesian_momentum_fixture=fixture,independent_viscosity_row_fixture=units,
                global_admissible_stress_lift_constructed=False,whole_outer_cone_certified=False,
                independently_bounded_global_flat_remainder=False,full_background_NS_validation=False,
                physical_energy_integral_certified=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual physical Gamma exterior PASS: source-bound map, full Cartesian momentum, viscosity units and regional exact zeros',flush=True)
    return result


if __name__=='__main__':run()
